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
from pymor.operators import LincombOperator
from pymor.parameters import CubicParameterSpace, Parametric
from pymor.reductors import reduce_generic_rb
from pymor.playground.algorithms import gram_schmidt_block_basis_extension, trivial_block_basis_extension
from pymor.playground.la import BlockVectorArray
from pymor.playground.reductors import GenericBlockRBReconstructor

from dune.pymor.core import wrap_module

from simdb.run import new_dataset, add_values, add_data, add_logfile

from OS2014_estimators import ReducedEstimator, reduce_with_estimator

import linearellipticexampleOS2014 as dune_module
DuneMatrixOperator = dune_module.Dune.Pymor.Operators.MatrixBasedDefault__DuneStuffLAEigenRowMajorSparseMatrix__lt___double___gt___DuneStuffLAEigenDenseVector__lt___double___gt__


config = {'dune_partitioning': '[8 8 1]',
          'dune_num_refinements': 2,
          'dune_oversampling_layers': 12,
          'dune_products': ['elliptic', 'penalty'],
          'dune_log_info_level': -1,
          'dune_log_debug_level': -1,
          'dune_log_enable_warnings': True,
          'dune_linear_solver_options': {'type': 'bicgstab.ilut', 'precision': '1e-14'},
          'dune_example': 'LocalThermalblockExample',
          'parameter_range': (0.1, 1.0),
          'compute_some_solution_norms': True,
          'initialize_basis_with': 1,
          'num_training_samples': 2,
          'greedy_max_extensions': 2,
          'extension_product': ('elliptic', 'penalty'),
          'greedy_error_norm': 'elliptic',
          'mu_hat_value': 1,
          'mu_bar_value': 1,
          'greedy_target_error': 1e-10,
          'greedy_use_estimator': True,
          'estimator_compute': 'model_reduction_error',
          'estimator_return': 'model_reduction_error',
          'num_test_samples': 11,
          'local_indicators': 'model_reduction_error',
          'marking_strategy': 'doerfler',
          'doerfler_marking_theta': 0.3,
          'local_boundary_values': 'dirichlet'}
DATASET_ID = config['dune_example'] + '_online_enrichment_test'

          # 'estimator_compute': ('discretization_error',
          #                       'full_error',
          #                       'model_reduction_error',
          #                       'eta_nc_red',
          #                       'eta_r_red',
          #                       'eta_df_red',
          #                       'eta_red'),


pymor.core.logger.MAX_HIERACHY_LEVEL = 2
getLogger('pymor.WrappedDiscretization').setLevel('WARN')
getLogger('pymor.algorithms').setLevel('INFO')
getLogger('dune.pymor.discretizations').setLevel('WARN')


def prepare(cfg):
    logger = getLogger('.OS2014.prepare')
    logger.setLevel('INFO')

    logger.info('Initializing DUNE module ({}):'.format(cfg['dune_example']))
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
    discretization = wrapper[example.discretization()]
    logger.info('  grid has {} subdomain{}'.format(discretization.num_subdomains,
                                                   '' if discretization.num_subdomains == 1 else 's'))
    logger.info('  discretization has {} DoFs'.format(discretization.solution_space.dim))
    discretization = discretization.with_(
            parameter_space=CubicParameterSpace(discretization.parameter_type,
                                                cfg['parameter_range'][0],
                                                cfg['parameter_range'][1]))
    logger.info('  parameter type is {}'.format(discretization.parameter_type))

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

    logger.info('Preparing products and norms ...')
    mu_hat_dune = wrapper.DuneParameter('mu', cfg['mu_hat_value'])
    mu_bar_dune = wrapper.DuneParameter('mu', cfg['mu_bar_value'])
    all_global_products = {}
    all_local_products = [{} for ss in np.arange(discretization.num_subdomains)]
    # get all products from the discretization and assemble the parametric ones
    for nm, pr in discretization.products.items():
        if not pr.parametric:
            all_global_products[nm] = pr
        else:
            all_global_products[nm] = pr.assemble(mu=wrapper[mu_bar_dune])
        for ss in np.arange(discretization.num_subdomains):
            pr = discretization.local_product(ss, nm)
            if pr.parametric:
                all_local_products[ss][nm] = pr.assemble(mu=wrapper[mu_bar_dune])
            else:
                all_local_products[ss][nm] = pr
    # combine products (if necessarry)
    global_product, _ = create_product(all_global_products, cfg['extension_product'])
    local_products, _ = map(list,
            zip(*[create_product(all_local_pr, cfg['extension_product']) for all_local_pr in all_local_products]))
    # assert all(nm == local_products_name[0] for nm in local_products_name)
    # local_products_name = local_products_name[0]
    norm, _ = create_product(all_global_products, cfg['greedy_error_norm'])
    norm = induced_norm(norm) if norm is not None else None
    return {'example': example,
            'discretization': discretization,
            'local_products': local_products,
            'norm': norm,
            'mu_bar_dune': mu_bar_dune,
            'mu_hat_dune': mu_hat_dune,
            'wrapper': wrapper}


def offline_phase(cfg, data):
    logger = getLogger('.OS2014.offline_phase')
    logger.setLevel('INFO')

    discretization = data['discretization']
    example        = data['example']
    local_products = data['local_products']
    norm           = data['norm']
    mu_bar_dune    = data['mu_bar_dune']
    mu_hat_dune    = data['mu_hat_dune']
    wrapper        = data['wrapper']

    # training_samples = list(discretization.parameter_space.sample_randomly(cfg['num_training_samples']))
    training_samples = list(discretization.parameter_space.sample_randomly(1))
    training_samples[0]['mu'][0] = 1.0
    add_values(training_samples=training_samples)
    if cfg['compute_some_solution_norms']:
        logger.info('Computing solution norms:')
        if norm is not None:
            solution_norms = [norm(discretization.solve(mu)) for mu in training_samples]
        else:
            solution_norms = [discretization.solve(mu).l2_norm()[0] for mu in training_samples]
        logger.info('  range:              [{}, {}]'.format(np.amin(solution_norms), np.amax(solution_norms)))
        logger.info('  median:              {}'.format(np.median(solution_norms)))
        logger.info('  mean:                {}'.format(np.mean(solution_norms)))
        logger.info('  standard deviation:  {}'.format(np.std(solution_norms)))
        add_values(solution_norms=solution_norms)
        logger.info('')

    extension_algorithm=partial(gram_schmidt_block_basis_extension, product=local_products)
    initial_basis = discretization.functionals['rhs'].source.empty()._blocks
    if cfg['initialize_basis_with'] >= 0:
        logger.info('Initializing local bases of up to order {} ...'.format(cfg['initialize_basis_with']))
        one = wrapper.vector_array(wrapper[example.project('1')])
        one = BlockVectorArray([discretization.localize_vector(one, ss)
                                for ss in np.arange(discretization.num_subdomains)])
        initial_basis, _ = extension_algorithm(initial_basis, one)
    if cfg['initialize_basis_with'] >= 1:
        xx = wrapper.vector_array(wrapper[example.project('x[0]')])
        xx = BlockVectorArray([discretization.localize_vector(xx, ss)
                               for ss in np.arange(discretization.num_subdomains)])
        initial_basis, _ = extension_algorithm(initial_basis, xx)
        yy = wrapper.vector_array(wrapper[example.project('x[1]')])
        yy = BlockVectorArray([discretization.localize_vector(yy, ss)
                               for ss in np.arange(discretization.num_subdomains)])
        initial_basis, _ = extension_algorithm(initial_basis, yy)
        xy = wrapper.vector_array(wrapper[example.project('x[0]*x[1]')])
        xy = BlockVectorArray([discretization.localize_vector(xy, ss)
                               for ss in np.arange(discretization.num_subdomains)])
        initial_basis, _ = extension_algorithm(initial_basis, xy)
    if cfg['initialize_basis_with'] >= 2:
        logger.warn('Ignoring initialize_basis_with higher than 1: {}'.format(cfg['initialize_basis_with']))
    # logger.info('  range (local basis size): [{}, {}]'.format(np.min([len(bb) for bb in initial_basis]),
    #                                                           np.max([len(bb) for bb in initial_basis])))

    reduced_estimator = ReducedEstimator(
            discretization, example, wrapper, mu_hat_dune, mu_bar_dune, norm, cfg['estimator_compute'], cfg['estimator_return'])
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
               offline_estimator_data=reduced_estimator.data)
    rd, rc = greedy_data['reduced_discretization'], greedy_data['reconstructor']
    basis, basis_mus = greedy_data['basis'], greedy_data['max_err_mus']
    return {'basis': basis,
            'basis_mus': basis_mus,
            'rc': rc,
            'rd': rd}


def online_phase(cfg, detailed_data, offline_data):
    logger = getLogger('.OS2014.offline_phase')
    logger.setLevel('INFO')

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

    def compute_oversampled_correction(oversampled_discretization, current_solution, mu, wrapper):
        mu_dune = wrapper.dune_parameter(mu)
        lhs = oversampled_discretization.operator._impl.freeze_parameter(mu_dune).container()
        rhs = oversampled_discretization.rhs._impl.freeze_parameter(mu_dune).as_vector()
        tmp = rhs.copy()
        lhs.mv(current_solution, tmp)
        rhs.isub(tmp)
        lhs_op = DuneMatrixOperator(lhs)
        a = b
        lhs_op.apply_inverse(rhs, tmp)
        return wrapper.vector_array(wrapper[tmp])

    discretization = detailed_data['discretization']
    example        = detailed_data['example']
    local_products = detailed_data['local_products']
    mu_bar_dune    = detailed_data['mu_bar_dune']
    mu_hat_dune    = detailed_data['mu_hat_dune']
    norm           = detailed_data['norm']
    wrapper        = detailed_data['wrapper']

    basis     = offline_data['basis']
    basis_mus = offline_data['basis_mus']
    rd        = offline_data['rd']
    rc        = offline_data['rc']

    reduced_estimator = ReducedEstimator(
            discretization, example, wrapper, mu_hat_dune, mu_bar_dune, norm, cfg['estimator_compute'], cfg['estimator_return'])
    reduced_estimator.extension_step += 1
    reduced_estimator.rc = rc
    num_test_samples = cfg['num_test_samples']
    target_error = cfg['greedy_target_error']

    logger.info('Started online phase on {} samples'.format(num_test_samples))
    # test_samples = list(discretization.parameter_space.sample_randomly(num_test_samples))
    test_samples = list(discretization.parameter_space.sample_randomly(1))
    test_samples[0]['mu'][0] = 0.1
    for mu in test_samples:
        mu_dune = wrapper.dune_parameter(mu)
        mu_in_basis = mu in basis_mus
        logger.info('Solving for {} ...'.format(mu))
        U_red = rd.solve(mu)
        logger.info('  Estimating (is {}in the basis):'.format('already ' if mu_in_basis else 'not '))
        error = reduced_estimator.estimate(U_red, mu, discretization)
        logger.info('              {}'.format(error))
        if error > target_error:
            if mu_in_basis:
                logger.error('Error ({}) is larger than target_error ({}), ' \
                             + 'but ({}) is already in the basis: aborting (this should not happen)!'.format(
                                 error, target_error, mu))
                break
            else:
                intermediate_basis = [bb.copy() for bb in basis]
                if cfg['local_indicators'] == 'model_reduction_error':
                    U_h = discretization.solve(mu)
                    assert len(U_h) == 1
                    U_h_dune = discretization.globalize_vectors(U_h)._list[0]._impl
                    discretization._impl.visualize(U_h_dune, 'detailed_solution_global', 'U_h')
                failure = False
                logger.info('Error ({}) is too large, starting intermediate offline phase:'.format(error))
                while error > target_error:
                    logger.info('  Estimating local error contributions ...')
                    U_red_h = rc.reconstruct(U_red)
                    assert len(U_red_h) == 1
                    U_red_global = discretization.globalize_vectors(U_red_h)
                    U_red_dune = U_red_global._list[0]._impl
                    discretization._impl.visualize(U_red_dune, 'reduced_solution_global', 'U_red')
                    # compute difference | U_h - U_red |
                    difference = U_h - U_red_h
                    for block in difference._blocks:
                        vec = block._list[0]._impl
                        for ii in np.arange(block.dim):
                            vec.set_entry(ii, np.abs(vec.get_entry(ii)))
                    difference_global = discretization.globalize_vectors(difference)
                    difference_dune = difference_global._list[0]._impl
                    discretization._impl.visualize(difference_dune, 'model_reduction_error_global', 'model reduction error')
                    # compute local error indicators
                    if cfg['local_indicators'] == 'model_reduction_error':
                        local_indicators = [None for ss in np.arange(discretization.num_subdomains)]
                        for subdomain in np.arange(discretization.num_subdomains):
                            local_norm = induced_norm(local_products[subdomain])
                            local_indicators[subdomain] = local_norm(difference.block(subdomain))
                    elif cfg['local_indicators'] == 'eta_red':
                        local_indicators = list(example.estimate_local(U_red_dune,
                                                                       'eta_OS2014_*',
                                                                       mu_hat_dune,
                                                                       mu_bar_dune,
                                                                       mu_dune))
                    else:
                        logger.error('Unknown local_indicators given: {}'.format(cfg['local_indicators']))
                        failure = True
                        break
                    # mark subdomains
                    if 'doerfler' in cfg['marking_strategy']:
                        marked_subdomains = doerfler_marking(local_indicators, cfg['doerfler_marking_theta'])
                    else:
                        logger.error('Unknown marking_strategy given: {}'.format(cfg['local_indicators']))
                        failure = True
                        break
                    logger.info('  {} subdomains marked, computing local solutions ...'.format(len(marked_subdomains)))
                    # compute updated local solution
                    local_solutions = [None for ss in np.arange(discretization.num_subdomains)]
                    for subdomain in marked_subdomains:
                        local_boundary_values = cfg['local_boundary_values']
                        if not (local_boundary_values == 'dirichlet' or local_boundary_values == 'neumann'):
                            logger.error('Unknown local_boundary_values given: {}'.format(local_boundary_values))
                            failure = True
                            break
                        oversampled_discretization = discretization.get_oversampled_discretization(
                                subdomain, local_boundary_values)
                        local_discretization = discretization.get_local_discretization(subdomain)
                        U_h_oversampled_dune   = example.project_global_to_oversampled(U_h_dune, subdomain)
                        U_h_oversampled = wrapper.vector_array(wrapper[U_h_oversampled_dune])
                        U_red_oversampled_dune = example.project_global_to_oversampled(U_red_dune, subdomain)
                        difference_oversampled_dune = example.project_global_to_oversampled(difference_dune, subdomain)
                        oversampled_discretization._impl.visualize(U_h_oversampled_dune,
                                                                   'detailed_solution_oversampled_{}'.format(subdomain),
                                                                   'U_h')
                        oversampled_discretization._impl.visualize(U_red_oversampled_dune,
                                                                   'reduced_solution_oversampled_{}'.format(subdomain),
                                                                   'U_red')
                        oversampled_discretization._impl.visualize(difference_oversampled_dune,
                                                                   'model_reduction_error_oversampled_{}'.format(subdomain),
                                                                   'model reduction error')
                        U_h_improved_oversampled_dune = example.solve_oversampled(
                                subdomain, local_boundary_values, U_red_oversampled_dune, mu_dune)
                        U_h_improved_oversampled = wrapper.vector_array(wrapper[U_h_improved_oversampled_dune])
                        oversampled_discretization._impl.visualize(U_h_improved_oversampled_dune,
                                                                   'improved_solution_oversampled_{}'.format(subdomain),
                                                                   'U_red')
                        difference_improved_oversampled = U_h_oversampled - U_h_improved_oversampled
                        difference_improved_oversampled_dune = difference_improved_oversampled._list[0]._impl
                        for ii in np.arange(difference_improved_oversampled_dune.dim()):
                            difference_improved_oversampled_dune.set_entry(ii, np.abs(difference_improved_oversampled_dune.get_entry(ii)))
                        oversampled_discretization._impl.visualize(difference_improved_oversampled_dune,
                                                                   'improved_model_reduction_error_oversampled_{}'.format(subdomain),
                                                                   'model reduction error')
                        U_red_local_dune = example.project_oversampled_to_local(U_red_oversampled_dune, subdomain)
                        U_red_local = wrapper.vector_array(wrapper[U_red_local_dune])
                        U_h_local_dune = example.project_global_to_local(U_h_dune, subdomain)
                        U_h_local = wrapper.vector_array(wrapper[U_h_local_dune])
                        difference_local = U_h_local - U_red_local
                        difference_local_dune = difference_local._list[0]._impl
                        for ii in np.arange(difference_local_dune.dim()):
                            difference_local_dune.set_entry(ii, np.abs(difference_local_dune.get_entry(ii)))
                        local_discretization._impl.visualize(difference_local_dune,
                                                             'model_reduction_error_local_{}'.format(subdomain),
                                                             'model reduction error')
                        local_norm = induced_norm(local_products[subdomain])
                        # print('|| U_red|_T - U_h|_T ||_T = {}'.format(local_norm(difference_local)[0]))
                        U_h_improved_local_dune = example.project_oversampled_to_local(
                                U_h_improved_oversampled_dune, subdomain)
                        U_h_improved_local = wrapper.vector_array(wrapper[U_h_improved_local_dune])
                        local_solutions[subdomain] = U_h_improved_local
                        difference_improved_local = U_h_local - U_h_improved_local
                        difference_improved_local_dune = difference_improved_local._list[0]._impl
                        for ii in np.arange(difference_improved_local_dune.dim()):
                            difference_improved_local_dune.set_entry(ii, np.abs(difference_improved_local_dune.get_entry(ii)))
                        local_discretization._impl.visualize(difference_improved_local_dune,
                                                             'improved_model_reduction_error_local_{}'.format(subdomain),
                                                             'model reduction error')
                        # print('|| U_h_improved - U_h|_T ||_T = {}'.format(local_norm(difference_improved_local)[0]))
                    # extend local bases
                    logger.info('  Enriching local bases on {} subdomain{}...'.format(
                        len(marked_subdomains), '' if len(marked_subdomains) == 1 else 's'))
                    try:
                        extended_bases, _ = gram_schmidt_block_basis_extension(
                                [basis[ss] for ss in marked_subdomains],
                                BlockVectorArray([local_solutions[ss] for ss in marked_subdomains], copy=False),
                                product=[local_products[ss] for ss in marked_subdomains])
                    except ExtensionError:
                        logger.error('Enrichment failed on all subdomains, aborting!')
                        failure = True
                        break
                    assert len(extended_bases) == len(marked_subdomains)
                    for ii, subdomain in enumerate(marked_subdomains):
                        intermediate_basis[subdomain] = extended_bases[ii]
                    logger.info('  Reducing ...')
                    rd, _, _ = reduce_generic_rb(discretization, intermediate_basis)
                    rc = GenericBlockRBReconstructor(intermediate_basis)
                    reduced_estimator.rc = rc
                    reduced_estimator.extension_step += 1
                    U_red = rd.solve(mu)
                    logger.info('  Estimating:')
                    new_error = reduced_estimator.estimate(U_red, mu, discretization)
                    logger.info('              {}'.format(new_error))
                    if new_error > error:
                        logger.error('Error increased, aborting!')
                        failure = True
                        break
                    error = new_error
                    if failure:
                        break
                if failure:
                    break
                logger.info('  Error ({}) is below the target error, continuing ...'.format(error))


if __name__ == '__main__':

    logfile = NamedTemporaryFile(delete=False).name
    pymor.core.logger.FILENAME = logfile
    new_dataset(DATASET_ID, **config)

    detailed_data = prepare(config)
    print('')
    offline_data  = offline_phase(config, detailed_data)
    print('')
    _             = online_phase(config, detailed_data, offline_data)

    # this should be the last action (to really capture all logs)
    add_logfile(logfile)

