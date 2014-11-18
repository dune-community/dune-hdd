#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from tempfile import NamedTemporaryFile
from functools import partial
import numpy as np

from pymor.algorithms import greedy
import pymor.core
from pymor.core import getLogger
from pymor.core.exceptions import ExtensionError
from pymor.la import induced_norm
from pymor.parameters import CubicParameterSpace, Parametric
from pymor.operators import LincombOperator
from pymor.reductors import reduce_generic_rb
from pymor.playground.la import BlockVectorArray
from pymor.playground.algorithms import gram_schmidt_block_basis_extension, trivial_block_basis_extension
from pymor.playground.reductors import GenericBlockRBReconstructor

from dune.pymor.core import wrap_module

from simdb.run import new_dataset, add_values, add_data, add_logfile

from OS2014_estimators import ReducedEstimator, reduce_with_estimator
import linearellipticexampleOS2014 as dune_module


dune_config = {'dune_partitioning': '[3 3 1]',
               'dune_num_refinements': 2,
               'dune_oversampling_layers': 21,
               'dune_products': ['l2', 'h1_semi', 'elliptic'],
               'dune_log_info_level': -1,
               'dune_log_debug_level': -1,
               'dune_log_enable_warnings': True,
               'dune_linear_solver_options': {'type': 'bicgstab.ilut', 'precision': '1e-14'},
               'dune_example': 'OS2014Example'}
config = {'num_training_samples': 10,
          'greedy_max_extensions': 11,
          'mu_hat_value': 0.1,
          'mu_bar_value': 0.1,
          'greedy_use_estimator': True,
          'greedy_target_error': 1e-14,
          'initialize_with_one': True,
          'estimator_return': 'eta_red',
          'num_test_samples': 20,
          'online_target_error': 0.275,
          'doerfler_marking_theta': 0.3}
DATASET_ID = dune_config['dune_example'] + '_online_enrichment_test'


pymor.core.logger.MAX_HIERACHY_LEVEL = 2
getLogger('pymor.WrappedDiscretization').setLevel('WARN')
getLogger('pymor.algorithms').setLevel('INFO')
getLogger('dune.pymor.discretizations').setLevel('WARN')


def get_logger():
    tmp_filename = NamedTemporaryFile(delete=False).name
    pymor.core.logger.FILENAME = tmp_filename
    logger = getLogger('.OS2014.main')
    logger.setLevel('INFO')
    return logger, tmp_filename


def init_dune(cfg):
    Example = dune_module.__dict__[cfg['dune_example']]
    example = Example(partitioning=cfg['dune_partitioning'],
                      num_refinements=cfg['dune_num_refinements'],
                      oversampling_layers=cfg['dune_oversampling_layers'],
                      products=cfg['dune_products'],
                      info_log_levels=cfg['dune_log_info_level'],
                      debug_log_levels=cfg['dune_log_debug_level'],
                      enable_warnings=cfg['dune_log_enable_warnings'],
                      enable_colors=True,
                      info_color='blue')
    _, wrapper = wrap_module(dune_module)
    return example, wrapper


def create_product(products, product_type):

    if product_type is None:
        return None, 'None'

    if not isinstance(product_type, tuple):
        return products[product_type], product_type
    else:
        prods = [products[tt] for tt in product_type]
        product_name = product_type[0]
        for ii in np.arange(1, len(product_type)):
            product_name += '_plus_' + product_type[ii]
        return LincombOperator(operators=prods, coefficients=[1 for pp in prods]), product_name


def create_discretization(example, wrapper, cfg):

    return wrapper[example.discretization()]
    # discretization = example.discretization()
    # linear_solver_options = cfg['dune_linear_solver_options']
    # dune_linear_solver_options = discretization.solver_options(linear_solver_options['type'])
    # for kk, vv in linear_solver_options.items():
    #     dune_linear_solver_options.set(kk, vv, True)
    # return wrapper[discretization].with_(solver_options=dune_linear_solver_options)


def compute_discretization_error(example, wrapper, discretization, norm_type, mu=None, mu_norm=None):
    U_h = discretization.solve(mu)
    assert len(U_h) == 1
    U_h_dune = U_h._list[0]._impl
    if mu_norm is not None:
        return example.compute_error(U_h_dune,
                                     norm_type if norm_type in dune_config['dune_products'] else 'l2',
                                     wrapper.dune_parameter(mu),
                                     mu_norm)
    else:
        return example.compute_error(U_h_dune,
                                     norm_type if norm_type in dune_config['dune_products'] else 'l2',
                                     wrapper.dune_parameter(mu))


def doerfler_marking(indicators, theta):
    assert 0.0 < theta <= 1.0
    indices = list(range(len(indicators)))
    # indicators = [ii**2 for ii in indicators]
    indicators, indices = [list(x) for x in zip(*sorted(zip(indicators, indices),
                                                        key=lambda pair: pair[0],
                                                        reverse=True))]
    total = np.sum(indicators)
    sums = np.array([np.sum(indicators[:ii+1]) for ii in np.arange(len(indicators))])
    where = sums > theta*total
    if np.any(where):
        return indices[:np.argmax(where)+1]
    else:
        return indices


def run_experiment(example, wrapper, discretization, cfg, product, norm):

    logger, logfile = get_logger()
    logger.info('computing solution norms:')

    mu_hat_dune = wrapper.DuneParameter('mu', cfg['mu_hat_value'])
    mu_bar_dune = wrapper.DuneParameter('mu', cfg['mu_bar_value'])
    mu_hat = wrapper[mu_hat_dune]
    mu_bar = wrapper[mu_bar_dune]
    all_global_products = {}
    all_local_products = [{} for ss in np.arange(discretization.num_subdomains)]
    for nm, pr in discretization.products.items():
        if not pr.parametric:
            all_global_products[nm] = pr
        else:
            all_global_products[nm] = pr.assemble(mu=wrapper[mu_bar_dune])
        for ss in np.arange(discretization.num_subdomains):
            pr = discretization.local_product(ss, nm)
            if pr.parametric:
                all_local_products[ss][nm] = pr.assemble(mu=mu_bar)
            else:
                all_local_products[ss][nm] = pr

    mu_norm = None
    if norm is not None and norm in dune_config['dune_products'] and discretization.products[norm].parametric:
        mu_norm = mu_bar_dune

    global_product, global_product_name = create_product(all_global_products, product)
    local_products, local_products_name = map(list,
                                              zip(*[create_product(all_local_pr, product)
                                                    for all_local_pr in all_local_products]))
    assert all(nm == local_products_name[0] for nm in local_products_name)
    local_products_name = local_products_name[0]
    norm, norm_name = create_product(all_global_products, norm)
    norm = induced_norm(norm) if norm is not None else None
    cfg['extension_product'] = local_products_name
    cfg['greedy_error_norm'] = norm_name

    new_dataset(DATASET_ID, **cfg)

    training_samples = list(discretization.parameter_space.sample_randomly(cfg['num_training_samples']))
    if norm is not None:
        solution_norms = [norm(discretization.solve(mu)) for mu in training_samples]
    else:
        solution_norms = [discretization.solve(mu).l2_norm()[0] for mu in training_samples]
    logger.info('  range:              [{}, {}]'.format(np.amin(solution_norms), np.amax(solution_norms)))
    logger.info('  median:              {}'.format(np.median(solution_norms)))
    logger.info('  mean:                {}'.format(np.mean(solution_norms)))
    logger.info('  standard deviation:  {}'.format(np.std(solution_norms)))
    add_values(training_samples=training_samples,
               solution_norms=solution_norms)
    logger.info('')

    # logger.info('computing some discretization errors ...')
    # less_samples = list(CubicParameterSpace(discretization.parameter_type, 0.1, 1.0).sample_randomly(10))
    # discretization_error_computer = partial(compute_discretization_error,
    #                                         example=example,
    #                                         wrapper=wrapper,
    #                                         discretization=discretization,
    #                                         norm_type=norm_name,
    #                                         mu_norm=mu_norm)
    # discretization_errors = [discretization_error_computer(mu=mu) for mu in less_samples]
    # logger.info('  range:              [{}, {}]'.format(np.amin(discretization_errors), np.amax(discretization_errors)))
    # logger.info('  median:              {}'.format(np.median(discretization_errors)))
    # logger.info('  mean:                {}'.format(np.mean(discretization_errors)))
    # logger.info('  standard deviation:  {}'.format(np.std(discretization_errors)))
    # add_values(less_samples=less_samples,
    #            discretization_errors=discretization_errors)
    # logger.info('')

    extension_algorithm=partial(gram_schmidt_block_basis_extension, product=local_products)
    initial_basis = discretization.functionals['rhs'].source.empty()._blocks
    if cfg['initialize_with_one']:
        one = wrapper.vector_array(wrapper[example.project('1')])
        one = BlockVectorArray([discretization.localize_vector(one, ss)
                                for ss in np.arange(discretization.num_subdomains)])
        initial_basis, _ = extension_algorithm(initial_basis, one)

    reduced_estimator = ReducedEstimator(discretization, example, wrapper, mu_hat_dune, mu_bar_dune, norm)
    greedy_data = greedy(discretization,
                         partial(reduce_with_estimator,
                                 reduced_estimator=reduced_estimator),
                         training_samples,
                         initial_basis=initial_basis,
                         use_estimator=cfg['greedy_use_estimator'],
                         error_norm=norm,
                         extension_algorithm=extension_algorithm,
                         max_extensions=cfg['greedy_max_extensions'],
                         target_error=cfg['greedy_target_error'])
    add_values(time=greedy_data['time'],
               max_err=greedy_data['max_err'],
               max_err_mus=greedy_data['max_err_mus'],
               extensions=greedy_data['extensions'],
               max_err_mu=greedy_data['max_err_mu'],
               max_errs=greedy_data['max_errs'],
               estimator_data=reduced_estimator.data)
    rd, rc = greedy_data['reduced_discretization'], greedy_data['reconstructor']
    basis, basis_mus = greedy_data['basis'], greedy_data['max_err_mus']

    # perform the 'online phase'
    num_test_samples = cfg['num_test_samples']
    target_error = cfg['online_target_error']
    print('')
    logger.info('Starting online phase for {} test parameters to reach a target error of {}:'.format(num_test_samples,
                                                                                                     target_error))
    test_samples = list(discretization.parameter_space.sample_randomly(num_test_samples))
    for mu in test_samples:
        mu_dune = wrapper.dune_parameter(mu)
        mu_in_basis = mu in basis_mus
        logger.info('Solving for {} (is {}in the basis) ...'.format(mu, 'already ' if mu_in_basis else 'not '))
        U_red = rd.solve(mu)
        logger.info('  Estimating ...')
        error = reduced_estimator.estimate(U_red, mu, discretization)
        if error > target_error:
            if mu_in_basis:
                logger.error('Error ({}) is larger than target_error ({}), ' \
                             + 'but ({}) is already in the basis: aborting!'.format(error, target_error, mu))
                break
            else:
                failure = False
                logger.info('Error ({}) is too large, starting intermediate offline phase:'.format(error))
                while error > target_error:
                    logger.info('  Estimating local error contributions ...')
                    U_red_h = rc.reconstruct(U_red)
                    assert len(U_red_h) == 1
                    U_red_global = discretization.globalize_vectors(U_red_h)
                    U_red_dune = U_red_global._list[0]._impl
                    local_indicators = list(example.estimate_local(U_red_dune,
                                                                   'eta_OS2014_*',
                                                                   mu_hat_dune,
                                                                   mu_bar_dune,
                                                                   mu_dune))
                    # marked_subdomains = doerfler_marking(local_indicators, cfg['doerfler_marking_theta'])
                    marked_subdomains = list(np.arange(discretization.num_subdomains))
                    logger.info('  Solving on {} subdomain{} ...'.format(len(marked_subdomains),
                                                                         '' if len(marked_subdomains) == 1 else 's'))
                    local_solutions = [None for ii in np.arange(len(marked_subdomains))]
                    for ii, subdomain in enumerate(marked_subdomains):
                        local_solutions[ii] = wrapper.vector_array(wrapper[example.solve_for_local_correction(
                                [U._list[0]._impl for U in U_red_h._blocks],
                                subdomain,
                                mu_dune)])
                    # local_solutions = discretization.solve(mu)._blocks
                    # local_solutions = [local_solutions[ss] for ss in marked_subdomains]
                    logger.info('  Enriching local bases on {} subdomain{}...'.format(
                        len(marked_subdomains), '' if len(marked_subdomains) == 1 else 's'))
                    try:
                        extended_bases, _ = gram_schmidt_block_basis_extension(
                                [basis[ss] for ss in marked_subdomains],
                                BlockVectorArray(local_solutions, copy=False),
                                product=[local_products[ss] for ss in marked_subdomains])
                    except ExtensionError:
                        logger.error('Enrichment failed on all subdomains, aborting!')
                        failure = True
                        break
                    assert len(extended_bases) == len(marked_subdomains)
                    for ii, subdomain in enumerate(marked_subdomains):
                        basis[subdomain] = extended_bases[ii]
                    logger.info('  Reducing ...')
                    rd, _, _ = reduce_generic_rb(discretization, basis)
                    rc = GenericBlockRBReconstructor(basis)
                    reduced_estimator.rc = rc
                    reduced_estimator.extension_step += 1
                    U_red = rd.solve(mu)
                    logger.info('  Estimating ...')
                    new_error = example.estimate(discretization.globalize_vectors(rc.reconstruct(U_red))._list[0]._impl, 'eta_OS2014_*', mu_hat_dune, mu_bar_dune, mu_dune)
                    # new_error = reduced_estimator.estimate(U_red, mu, discretization)
                    if new_error > error:
                        logger.error('Error increased, aborting!')
                        failure = True
                        break
                    error = new_error
                if failure:
                    break
                logger.info('  Error ({}) is below the target error, continuing ...'.format(error))

    # this should be the last action (to really capture all logs)
    add_logfile(logfile)


if __name__ == '__main__':

    print('initializing dune module... ')
    example, wrapper = init_dune(dune_config)
    discretization = create_discretization(example, wrapper, dune_config)
    discretization = discretization.with_(
                    parameter_space=CubicParameterSpace(discretization.parameter_type, 0.1, 1.0))
    print('the discretization has {} DoFs and {} subdomains.'.format(discretization.solution_space.dim,
                                                                     discretization.num_subdomains))
    print('the parameter type is {}.'.format(discretization.parameter_type))
    print('')

    product = ('l2', 'h1_semi')
    norm = 'elliptic'
    cfg = config.copy()
    for kk, vv in dune_config.items():
        assert not cfg.has_key(kk)
        cfg[kk] = vv

    print('running experiment with product \'{}\' and norm \'{}\':'.format(product, norm))
    run_experiment(example, wrapper, discretization, cfg, product, norm)
    print('')

