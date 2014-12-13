#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from functools import partial
import numpy as np

from pymor.algorithms import greedy
from pymor.core import getLogger
from pymor.playground.algorithms import gram_schmidt_block_basis_extension
from pymor.playground.la import BlockVectorArray

from simdb.run import add_values

from OS2014_estimators import ReducedEstimator, reduce_with_estimator


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

    if cfg['training_sampling_strategy'] == 'random':
        training_samples = list(discretization.parameter_space.sample_randomly(cfg['num_training_samples']))
    elif cfg['training_sampling_strategy'] == 'uniform':
        training_samples = list(discretization.parameter_space.sample_uniformly(cfg['num_training_samples']))
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

