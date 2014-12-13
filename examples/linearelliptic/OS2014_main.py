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

from OS2014_estimators import DetailedEstimator, ReducedEstimator, reduce_with_estimator

import linearellipticexampleOS2014 as dune_module


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
          'compute_some_solution_norms': False,
          'initialize_basis_with': 1,
          'num_training_samples': 0,
          'greedy_max_extensions': 0,
          'extension_product': ('elliptic', 'penalty'),
          'greedy_error_norm': 'elliptic',
          'mu_hat_value': 1,
          'mu_bar_value': 1,
          'greedy_target_error': 1e-10,
          'greedy_use_estimator': True,
          'estimator_compute': 'model_reduction_error',
          'estimator_return': 'model_reduction_error',
          'num_test_samples': 10,
          'estimate_some_errors': False,
          'local_indicators': 'model_reduction_error',
          'marking_strategy': 'doerfler_and_age',
          'doerfler_marking_theta': 0.75,
          'local_boundary_values': 'dirichlet',
          'online_target_error': 1e-4,
          'online_max_extensions': 20}
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


class ConfigurationError(Exception):
    """Raised when an invalid configuration is given."""


class EnrichmentError(Exception):
    """Raised when the enrichment failed."""


def prepare(cfg):
    logger = getLogger('.OS2014.prepare')
    logger.setLevel('INFO')
    reference_needed = ((cfg['greedy_use_estimator'] or cfg['estimate_some_errors'])
                         and ('discretization_error' in cfg['estimator_compute']
                              or 'full_error' in cfg['estimator_compute']
                              or 'model_reduction_error' in cfg['estimator_compute'])
                        or cfg['local_indicators'] == 'model_reduction_error')

    logger.info('Initializing DUNE module ({}):'.format(cfg['dune_example']))
    Example = dune_module.__dict__[cfg['dune_example']]
    example = Example(partitioning=cfg['dune_partitioning'],
                      num_refinements=cfg['dune_num_refinements'],
                      oversampling_layers=cfg['dune_oversampling_layers'],
                      products=cfg['dune_products'],
                      with_reference=reference_needed,
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

    training_samples = list(discretization.parameter_space.sample_randomly(cfg['num_training_samples']))
    add_values(training_samples=training_samples)
    if cfg['compute_some_solution_norms'] and len(training_samples) > 0:
        logger.info('Computing solution norms:')
        if norm is not None:
            solution_norms = [norm(discretization.solve(mu)) for mu in training_samples]
        else:
            solution_norms = [discretization.solve(mu).l2_norm()[0] for mu in training_samples]
        logger.info('  range:              [{}, {}]'.format(np.amin(solution_norms), np.amax(solution_norms)))
        logger.info('  mean:                {}'.format(np.mean(solution_norms)))
        add_values(solution_norms=solution_norms)

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
    logger.info('')

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
               max_err_mus=greedy_data['max_err_mus'],
               extensions=greedy_data['extensions'],
               max_errs=greedy_data['max_errs'],
               offline_estimator_data=reduced_estimator.data)
    rd, rc = greedy_data['reduced_discretization'], greedy_data['reconstructor']
    basis, basis_mus = greedy_data['basis'], greedy_data['max_err_mus']

    print('')
    logger.info('Offline phase finished.')
    logger.info('Basis sizes range from {} to {}.'.format(np.min([len(bb) for bb in basis]),
                                                          np.max([len(bb) for bb in basis])))
    return {'basis': basis,
            'basis_mus': basis_mus,
            'rc': rc,
            'rd': rd}


def online_phase(cfg, detailed_data, offline_data):
    logger = getLogger('.OS2014.online_phase')
    logger.setLevel('INFO')

    def doerfler_marking(indicators, theta):
        assert 0.0 < theta <= 1.0
        indices = list(range(len(indicators)))
        indicators = [ii**2 for ii in indicators]
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
    target_error = cfg['online_target_error']

    logger.info('Started online phase for {} samples'.format(num_test_samples))
    test_samples = list(discretization.parameter_space.sample_randomly(num_test_samples))

    if cfg['estimate_some_errors'] and len(test_samples) > 0:
        logger.info('Estimating discretization errors:')
        detailed_estimator = DetailedEstimator(example, wrapper, mu_hat_dune, mu_bar_dune)
        estimates = [detailed_estimator.estimate(
                         discretization.globalize_vectors(discretization.solve(mu))._list[0]._impl,
                         wrapper.dune_parameter(mu)) for mu in test_samples]
        max_error = np.amax(estimates)
        logger.info('  range:              [{}, {}]'.format(np.amin(estimates), max_error))
        logger.info('  mean:                {}'.format(np.mean(estimates)))
        add_values(estimates=estimates)
        if max_error > cfg['online_target_error']:
            logger.warn('Given target error of {} is below the worst discretization error {}!'.format(
                cfg['online_target_error'], max_error))
        print('')

    failures = 0
    successes = 0
    age = np.ones(discretization.num_subdomains)
    for mu in test_samples:
        mu_dune = wrapper.dune_parameter(mu)
        mu_in_basis = mu in basis_mus
        logger.info('Solving for {} ...'.format(mu))
        U_red = rd.solve(mu)
        logger.info('Estimating (mu is {}in the basis) ...'.format('already ' if mu_in_basis else 'not '))
        error = reduced_estimator.estimate(U_red, mu, discretization)
        if error > target_error:
            if mu_in_basis:
                logger.error(('Error ({}) is larger than the target_error ({}), '
                             + 'but {} is already in the basis: aborting!').format(
                                 error, target_error, mu))
                logger.error('This usually means that the tolerances are poorly chosen!')
                failures += 1
                break
            else:
                try:
                    logger.info('Error ({}) is too large, starting local enrichment phase:'.format(error))
                    num_extensions = 0

                    intermediate_basis = [bb.copy() for bb in basis]
                    if cfg['local_indicators'] == 'model_reduction_error':
                        U_h = discretization.solve(mu)
                        assert len(U_h) == 1
                    while error > target_error and num_extensions < cfg['online_max_extensions']:
                        logger.info('  Estimating local error contributions ...')
                        U_red_h = rc.reconstruct(U_red)
                        assert len(U_red_h) == 1
                        U_red_global = discretization.globalize_vectors(U_red_h)
                        U_red_dune = U_red_global._list[0]._impl
                        # compute local error indicators
                        if cfg['local_indicators'] == 'model_reduction_error':
                            difference = U_h - U_red_h
                            local_indicators = [induced_norm(local_products[ss])(difference._blocks[ss])
                                                for ss in np.arange(discretization.num_subdomains)]
                        elif cfg['local_indicators'] == 'eta_red':
                            local_indicators = list(example.estimate_local(U_red_dune,
                                                                           'eta_OS2014_*',
                                                                           mu_hat_dune,
                                                                           mu_bar_dune,
                                                                           mu_dune))
                        else:
                            raise ConfigurationError('Unknown local_indicators given: {}'.format(cfg['local_indicators']))
                        # mark subdomains
                        if 'doerfler' in cfg['marking_strategy']:
                            marked_subdomains = set(doerfler_marking(local_indicators, cfg['doerfler_marking_theta']))
                        else:
                            raise ConfigurationError('Unknown marking_strategy given: {}'.format(cfg['local_indicators']))
                        if 'neighbours' in cfg['marking_strategy']:
                            for ss in list(marked_subdomains):
                                neighbours = (list(discretization._impl.neighbouring_subdomains(ss)))
                                for nn in neighbours:
                                    marked_subdomains.append(nn)
                            marked_subdomains = set(marked_subdomains)
                        if 'age' in cfg['marking_strategy']:
                            only_marked = len(marked_subdomains)
                            too_old = np.where(age > cfg['marking_max_age'])[0]
                            for ss in too_old:
                                marked_subdomains.add(ss)
                            logger.info('  {} subdomains marked ({} bc. of age), computing local solutions ...'.format(
                                len(marked_subdomains), len(marked_subdomains) - only_marked))
                        else:
                            logger.info('  {} subdomains marked, computing local solutions ...'.format(
                                len(marked_subdomains)))
                        for ss in np.arange(discretization.num_subdomains):
                            if ss in marked_subdomains:
                                age[ss] = 0
                            else:
                                age[ss] += 1
                        # compute updated local solution
                        local_solutions = [None for ss in np.arange(discretization.num_subdomains)]
                        for subdomain in marked_subdomains:
                            local_boundary_values = cfg['local_boundary_values']
                            if not (local_boundary_values == 'dirichlet' or local_boundary_values == 'neumann'):
                                raise ConfigurationError('Unknown local_boundary_values given: {}'.format(local_boundary_values))
                            oversampled_discretization = discretization.get_oversampled_discretization(
                                    subdomain, local_boundary_values)
                            local_discretization = discretization.get_local_discretization(subdomain)
                            U_red_oversampled_dune = example.project_global_to_oversampled(U_red_dune, subdomain)
                            U_h_improved_oversampled_dune = example.solve_oversampled(
                                    subdomain, local_boundary_values, U_red_oversampled_dune, mu_dune)
                            U_h_improved_local_dune = example.project_oversampled_to_local(
                                    U_h_improved_oversampled_dune, subdomain)
                            U_h_improved_local = wrapper.vector_array(wrapper[U_h_improved_local_dune])
                            local_solutions[subdomain] = U_h_improved_local
                        # extend local bases
                        logger.info('  Extending bases on {} subdomain{}...'.format(
                            len(marked_subdomains), '' if len(marked_subdomains) == 1 else 's'))
                        old_basis_size = sum([len(bb) for bb in intermediate_basis])
                        extended_bases, _ = gram_schmidt_block_basis_extension(
                                [intermediate_basis[ss] for ss in marked_subdomains],
                                [local_solutions[ss] for ss in marked_subdomains],
                                product=[local_products[ss] for ss in marked_subdomains])
                        assert len(extended_bases) == len(marked_subdomains)
                        for ii, subdomain in enumerate(marked_subdomains):
                            intermediate_basis[subdomain] = extended_bases[ii]
                        new_basis_size = sum([len(bb) for bb in intermediate_basis])
                        num_extensions += 1
                        logger.info('  Reducing ...')
                        rd, _, _ = reduce_generic_rb(discretization, intermediate_basis)
                        rc = GenericBlockRBReconstructor(intermediate_basis)
                        reduced_estimator.rc = rc
                        reduced_estimator.extension_step += 1
                        U_red = rd.solve(mu)
                        logger.info('  Estimating (total basis size: {})'.format(
                            sum(len(bb) for bb in intermediate_basis)))
                        new_error = reduced_estimator.estimate(U_red, mu, discretization)
                        order = np.log(new_error/error)/np.log(old_basis_size/new_basis_size)
                        logger.info('              {} (order: {})'.format(new_error, order))
                        if new_error > error:
                            logger.warn('Error increased (from {} to {}) after enrichment!'.format(error, new_error))
                        elif order < 1:
                            logger.warn('Error decreased only slightly (from {} to {}) after enrichment!'.format(error, new_error))
                        if num_extensions >= cfg['online_max_extensions']:
                            basis = intermediate_basis
                            raise EnrichmentError('Reached maximum number of {} extensions!'.format(
                                cfg['online_max_extensions']))
                        error = new_error
                    logger.info('  Error ({}) is below the target error, continuing ...'.format(error))
                    successes += 1
                    basis = intermediate_basis
                    logger.info('Basis sizes range from {} to {}.'.format(np.min([len(bb) for bb in basis]),
                                                                          np.max([len(bb) for bb in basis])))
                except EnrichmentError, ee:
                    logger.critical('Enrichment stopped because: {}'.format(ee))
                    logger.info('Basis sizes range from {} to {}.'.format(np.min([len(bb) for bb in basis]),
                                                                          np.max([len(bb) for bb in basis])))
                    logger.info('Continuing with the next parameter ...')
                    failures += 1
                print('')
        else:
            logger.info('Error ({}) is below the target error, continuing ...'.format(error))
            successes += 1
            print('')

    logger.info('Adaptive online phase finished.')
    if failures == 0 and len(test_samples) > 0:
        logger.info('  Target error could be reached for all {} parameters.'.format(len(test_samples)))
    elif successes == 0 and len(test_samples) > 0:
        logger.warn('  Target error could not be reached for any of the {} parameters!'.format(len(test_samples)))
    else:
        if successes > 0:
            logger.info('  Target error could be reached for {} out of {} parameters.'.format(successes, len(test_samples)))
        if failures > 0:
            logger.info('  Target error could not be reached for {} out of {} parameters.'.format(failures,
                                                                                              len(test_samples)))
    logger.info('Final Basis sizes range from {} to {}.'.format(np.min([len(bb) for bb in basis]),
                                                                np.max([len(bb) for bb in basis])))
    example.visualize_on_coarse_grid([len(bb) for bb in basis], 'final_basis_sizes', 'local_basis_size')


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

