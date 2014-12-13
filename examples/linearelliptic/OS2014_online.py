#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

import numpy as np

from pymor.core import getLogger
from pymor.la import induced_norm
from pymor.reductors import reduce_generic_rb
from pymor.playground.algorithms import gram_schmidt_block_basis_extension
from pymor.playground.reductors import GenericBlockRBReconstructor

from simdb.run import add_values

from OS2014_estimators import DetailedEstimator, ReducedEstimator


class ConfigurationError(Exception):
    """Raised when an invalid configuration is given."""


class EnrichmentError(Exception):
    """Raised when the enrichment failed."""


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
                logger.error(('The error ({}) is larger than the target_error ({}), '
                             + 'but {} is already in the basis: aborting!').format(
                                 error, target_error, mu))
                logger.error('This usually means that the tolerances are poorly chosen!')
                failures += 1
                print('')
            else:
                try:
                    logger.info('The error ({}) is too large, starting local enrichment phase:'.format(error))
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
                            logger.warn('The error has increased (from {} to {}) after enrichment!'.format(error, new_error))
                        elif order < 1:
                            logger.warn('The error has decreased only slightly (from {} to {}) after enrichment!'.format(error, new_error))
                        if num_extensions >= cfg['online_max_extensions'] and new_error > cfg['online_target_error']:
                            basis = intermediate_basis
                            raise EnrichmentError('Reached maximum number of {} extensions!'.format(
                                cfg['online_max_extensions']))
                        error = new_error
                    logger.info('  The error ({}) is below the target error, continuing ...'.format(error))
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
            logger.info('The error ({}) is below the target error, continuing ...'.format(error))
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
    final_basis_sizes = [len(bb) for bb in basis]
    add_values(final_basis_sizes=final_basis_sizes)
    example.visualize_on_coarse_grid(final_basis_sizes, 'final_basis_sizes', 'local_basis_size')

