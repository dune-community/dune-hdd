#!/usr/bin/env python
# This file is part of the dune-pymor project:
#   https://github.com/pyMor/dune-pymor
# Copyright Holders: Felix Albrecht, Stephan Rave
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

'''Dune LRBMS demo.

Usage:
  thermalblock_main.py CONFIGFILE

Arguments:
  CONFIGFILE File that can be understood by pythons ConfigParser and by dunes ParameterTree
'''

from __future__ import absolute_import, division, print_function

config_defaults = {'training_set': 'random',
                   'num_training_samples': '100',
                   'reductor': 'generic',
                   'extension_algorithm': 'gram_schmidt',
                   'extension_algorithm_product': 'h1_semi',
                   'greedy_error_norm': 'h1_semi',
                   'use_estimator': 'False',
                   'max_rb_size': '100',
                   'target_error': '0.01',
                   'final_compression': 'False',
                   'compression_product': 'None',
                   'test_set': 'training',
                   'num_test_samples': '100',
                   'test_error_norm': 'h1_semi'}

import sys
import math as m
import time
from functools import partial
from itertools import izip
import numpy as np
from docopt import docopt
from scipy.sparse import coo_matrix
from scipy.sparse import bmat as sbmat
from numpy import bmat as nbmat
import ConfigParser

import linearellipticexamplethermalblock as dune_module
from dune.pymor.core import wrap_module

import pymor.core as core
core.logger.MAX_HIERACHY_LEVEL = 2
from pymor import defaults
from pymor.algorithms import greedy, gram_schmidt_basis_extension, pod_basis_extension, trivial_basis_extension
from pymor.playground.algorithms import greedy_lrbms
from pymor.core import cache
from pymor.core.exceptions import ConfigError
from pymor.discretizations import StationaryDiscretization
from pymor.la import NumpyVectorArray
from pymor.la.basic import induced_norm
from pymor.la.pod import pod
from pymor.la import induced_norm
from pymor.operators import NumpyMatrixOperator
from pymor.operators.basic import NumpyLincombMatrixOperator
from pymor.parameters import CubicParameterSpace
from pymor.reductors import reduce_generic_rb
from pymor.reductors.basic import GenericRBReconstructor, reduce_generic_rb
from pymor.reductors.linear import reduce_stationary_affine_linear

logger = core.getLogger('pymor.main.demo')
logger.setLevel('INFO')
core.getLogger('pymor.WrappedDiscretization').setLevel('WARN')
core.getLogger('pymor.algorithms').setLevel('INFO')
core.getLogger('dune.pymor.discretizations').setLevel('WARN')

def load_dune_module(settings_filename):

    logger.info('initializing dune module...')
    example = dune_module.ThermalblockExample()
    example.initialize([settings_filename])
    _, wrapper = wrap_module(dune_module)
    return example, wrapper


def perform_standard_rb(config, detailed_discretization, training_samples):

    # parse config
    reductor_id = config.get('pymor', 'reductor')
    if reductor_id == 'generic':
        reductor = reduce_generic_rb
    elif reductor_id == 'stationary_affine_linear':
        reductor_error_product = config.get('pymor', 'reductor_error_product')
        assert reductor_error_product == 'None'
        reductor_error_product = None
        reductor = partial(reduce_stationary_affine_linear, error_product=reductor_error_product)
    else:
        raise ConfigError('unknown \'pymor.reductor\' given: \'{}\''.format(reductor_id))

    # first the extension algorithm product, if needed
    extension_algorithm_id = config.get('pymor', 'extension_algorithm')
    if extension_algorithm_id in {'gram_schmidt', 'pod'}:
        extension_algorithm_product_id = config.get('pymor', 'extension_algorithm_product')
        if extension_algorithm_product_id == 'None':
            extension_algorithm_product = None
        else:
            extension_algorithm_product = detailed_discretization.products[extension_algorithm_product_id]
    # then the extension algorithm
    extension_algorithm_id = config.get('pymor', 'extension_algorithm')
    if extension_algorithm_id == 'gram_schmidt':
        extension_algorithm = partial(gram_schmidt_basis_extension, product=extension_algorithm_product)
        extension_algorithm_id += ' ({})'.format(extension_algorithm_product_id)
    elif extension_algorithm_id == 'pod':
        extension_algorithm = partial(pod_basis_extension, product=extension_algorithm_product)
        extension_algorithm_id += ' ({})'.format(extension_algorithm_product_id)
    elif extension_algorithm_id == 'trivial':
        extension_algorithm = trivial_basis_extension
    else:
        raise ConfigError('unknown \'pymor.extension_algorithm\' given: \'{}\''.format(extension_algorithm_id))

    greedy_error_norm_id = config.get('pymor', 'greedy_error_norm')
    assert greedy_error_norm_id in {'l2', 'h1_semi'}
    greedy_error_norm = induced_norm(detailed_discretization.products[greedy_error_norm_id])

    greedy_use_estimator = config.getboolean('pymor', 'use_estimator')
    greedy_max_rb_size = config.getint('pymor', 'max_rb_size')
    greedy_target_error = config.getfloat('pymor', 'target_error')

    # do the actual work
    greedy_data = greedy(detailed_discretization,
                         reductor,
                         training_samples,
                         initial_basis=detailed_discretization.functionals['rhs'].type_source.empty(
                                      dim=detailed_discretization.functionals['rhs'].dim_source),
                         use_estimator=greedy_use_estimator,
                         error_norm=greedy_error_norm,
                         extension_algorithm=extension_algorithm,
                         max_extensions=greedy_max_rb_size,
                         target_error=greedy_target_error)
    reduced_basis = greedy_data['basis']
    rb_size = len(reduced_basis)

    # perform final compression
    final_compression = config.getboolean('pymor', 'final_compression')
    if final_compression:
        t = time.time()
        logger.info('Applying final compression ...')

        # select product
        compression_product_id = config.get('pymor', 'compression_product')
        if compression_product_id == 'None':
            compression_product = None
        else:
            compression_product = detailed_discretization.products[compression_product_id]

        # do the actual work
        reduced_basis = pod(reduced_basis, product=compression_product)

        rd, rc, _ = reductor(detailed_discretization, reduced_basis)
        greedy_data['reduced_discretization'], greedy_data['reconstructor'] = rd, rc

        time_compression = time.time() - t
        compressed_rb_size = len(reduced_basis)

        #report
        report_string = '''
Greedy basis generation:
    used estimator:        {greedy_use_estimator}
    error norm:            {greedy_error_norm_id}
    extension method:      {extension_algorithm_id}
    prescribed basis size: {greedy_max_rb_size}
    prescribed error:      {greedy_target_error}
    actual basis size:     {rb_size}
    greedy time:           {greedy_data[time]}
    compression method:     pod ({compression_product_id})
    compressed basis size:  {compressed_rb_size}
    final compression time: {time_compression}
'''.format(**locals())
    else:
        #report
        report_string = '''
Greedy basis generation:
    used estimator:        {greedy_use_estimator}
    error norm:            {greedy_error_norm_id}
    extension method:      {extension_algorithm_id}
    prescribed basis size: {greedy_max_rb_size}
    prescribed error:      {greedy_target_error}
    actual basis size:     {rb_size}
    elapsed time:          {greedy_data[time]}
'''.format(**locals())

    return report_string, greedy_data


def perform_lrbms(config, multiscale_discretization, training_samples):

    num_subdomains = multiscale_discretization._impl.num_subdomains()

    # parse config
    # first the extension algorithm product, if needed
    extension_algorithm_id = config.get('pymor', 'extension_algorithm')
    if extension_algorithm_id in {'gram_schmidt', 'pod'}:
        extension_algorithm_product_id = config.get('pymor', 'extension_algorithm_product')
        if extension_algorithm_product_id == 'None':
            extension_algorithm_products = [None for ss in np.arange(num_subdomains)]
        else:
            extension_algorithm_products = [multiscale_discretization.local_product(ss, extension_algorithm_product_id)
                                            for ss in np.arange(num_subdomains)]
    # then the extension algorithm
    if extension_algorithm_id == 'gram_schmidt':
        extension_algorithm = [partial(gram_schmidt_basis_extension, product=extension_algorithm_products[ss])
                               for ss in np.arange(num_subdomains)]
        extension_algorithm_id += ' ({})'.format(extension_algorithm_product_id)
    elif extension_algorithm_id == 'pod':
        extension_algorithm = [partial(pod_basis_extension, product=extension_algorithm_products[ss])
                               for ss in np.arange(num_subdomains)]
        extension_algorithm_id += ' ({})'.format(extension_algorithm_product_id)
    elif extension_algorithm_id == 'trivial':
        extension_algorithm = [trivial_basis_extension for ss in np.arange(num_subdomains)]
    else:
        raise ConfigError('unknown \'pymor.extension_algorithm\' given:\'{}\''.format(extension_algorithm_id))

    greedy_error_norm_id = config.get('pymor', 'greedy_error_norm')
    if greedy_error_norm_id == 'None':
        greedy_error_norm = None
    else:
        greedy_error_norm = induced_norm(multiscale_discretization.products[greedy_error_norm_id])

    greedy_use_estimator = config.getboolean('pymor', 'use_estimator')
    assert greedy_use_estimator is False
    greedy_max_rb_size = config.getint('pymor', 'max_rb_size')
    greedy_target_error = config.getfloat('pymor', 'target_error')

    # do the actual work
    greedy_data = greedy_lrbms(multiscale_discretization,
                               reduce_generic_rb,
                               training_samples,
                               initial_basis=[multiscale_discretization.local_rhs(ss).type_source.empty(dim=multiscale_discretization.local_rhs(ss).dim_source)
                                             for ss in np.arange(num_subdomains)],
                               use_estimator=greedy_use_estimator,
                               error_norm=greedy_error_norm,
                               extension_algorithm=extension_algorithm,
                               max_extensions=greedy_max_rb_size,
                               target_error=greedy_target_error)

    reduced_basis = greedy_data['basis']
    rb_size = [len(local_data) for local_data in reduced_basis]

    # perform final compression
    final_compression = config.getboolean('pymor', 'final_compression')
    if final_compression:
        t = time.time()
        logger.info('Applying final compression ...')

        # select local product
        compression_product_id = config.get('pymor', 'compression_product')
        if compression_product_id == 'None':
            compression_products = [None for ss in np.arange(num_subdomains)]
        else:
            compression_products = [multiscale_discretization.local_product(ss, compression_product_id)
                                            for ss in np.arange(num_subdomains)]

        # do the actual work
        reduced_basis = [pod(reduced_basis[ss], product=compression_products[ss]) for ss in np.arange(num_subdomains)]

        rd, rc, _ = reduce_generic_rb(multiscale_discretization, reduced_basis)
        greedy_data['reduced_discretization'], greedy_data['reconstructor'] = rd, rc

        time_compression = time.time() - t
        compressed_rb_size = [len(local_data) for local_data in reduced_basis]

        #report
        report_string = '''
Greedy basis generation:
    used estimator:        {greedy_use_estimator}
    error norm:            {greedy_error_norm_id}
    extension method:      {extension_algorithm_id}
    prescribed basis size: {greedy_max_rb_size}
    prescribed error:      {greedy_target_error}
    actual basis size:     {rb_size}
    greedy time:           {greedy_data[time]}
    compression method:     pod ({compression_product_id})
    compressed basis size:  {compressed_rb_size}
    final compression time: {time_compression}
'''.format(**locals())
    else:
        #report
        report_string = '''
Greedy basis generation:
    used estimator:        {greedy_use_estimator}
    error norm:            {greedy_error_norm_id}
    extension method:      {extension_algorithm_id}
    prescribed basis size: {greedy_max_rb_size}
    prescribed error:      {greedy_target_error}
    actual basis size:     {rb_size}
    elapsed time:          {greedy_data[time]}
'''.format(**locals())

    return report_string, greedy_data


def test_quality(config, test_samples, detailed_discretization, greedy_data, strategy = 'stochastic'):

    # parse config
    test_error_norm_id = config.get('pymor', 'test_error_norm')
    if test_error_norm_id == 'None':
        test_error_norm = None
    else:
        test_error_norm = induced_norm(detailed_discretization.products[test_error_norm_id])

    # get reduced quantities
    reduced_discretization = greedy_data['reduced_discretization']
    reconstructor          = greedy_data['reconstructor']

    # run the test
    test_size = len(test_samples)
    tic = time.time()

    def compute_error(mu):
        detailed_solution = detailed_discretization.solve(mu)
        reduced_DoF_vector = reduced_discretization.solve(mu)
        reduced_solution = reconstructor.reconstruct(reduced_DoF_vector)
        if test_error_norm is None:
            return (detailed_solution - reduced_solution).l2_norm()
        else:
            return test_error_norm(detailed_solution - reduced_solution)

    toc = time.time()

    errors = np.array([compute_error(mu) for mu in test_samples])
    max_err = -1
    maximizing_mu = []
    for ii in np.arange(test_size):
        if errors[ii] > max_err:
            max_err = errors[ii]
            maximizing_mu = test_samples[ii]
    max_err = max_err[0]

    median = np.median(errors)
    mean = np.mean(errors)
    stddev = np.std(errors)
    min_err = np.amin(errors)
    range_err = max_err - min_err

    t_est = toc - tic

    # and report
    return '''
{strategy} error estimation:
    error norm: {test_error_norm_id}
    number of samples: {test_size}
    elapsed time:      {t_est}
    maximizing mu:     {maximizing_mu}
    error
      range:              [{min_err}, {max_err}]
      median:             {median}
      mean:               {mean}
      standard deviation: {stddev}
'''.format(**locals())


if __name__ == '__main__':
    # first of all, clear the cache
    cache.clear_caches()
    # parse arguments
    args = docopt(__doc__)
    config = ConfigParser.ConfigParser(config_defaults)
    try:
        config.readfp(open(args['CONFIGFILE']))
        assert config.has_section('pymor')
    except:
        raise ConfigError('CONFIGFILE has to be the name of an existing file that contains a [pymor] section')
    if config.has_section('pymor.defaults'):
        float_suffixes = ['_tol', '_threshold']
        boolean_suffixes = ['_find_duplicates', '_check', '_symmetrize', '_orthonormalize', '_raise_negative',
                            'compact_print']
        int_suffixes = ['_maxiter']
        for key, value in config.items('pymor.defaults'):
            if any([len(key) >= len(suffix) and key[-len(suffix):] == suffix for suffix in float_suffixes]):
                defaults.__setattr__(key, config.getfloat('pymor.defaults', key))
            elif any([len(key) >= len(suffix) and key[-len(suffix):] == suffix for suffix in boolean_suffixes]):
                defaults.__setattr__(key, config.getboolean('pymor.defaults', key))
            elif any([len(key) >= len(suffix) and key[-len(suffix):] == suffix for suffix in int_suffixes]):
                defaults.__setattr__(key, config.getint('pymor.defaults', key))

    # load module
    example, wrapper = load_dune_module(args['CONFIGFILE'])

    # create global cg discretization
    discretization = wrapper[example.discretization()]
    discretization = discretization.with_(
        parameter_space=CubicParameterSpace(discretization.parameter_type, 0.1, 10.0))
    logger.info('The discretization has {} DoFs.'.format(discretization.operator.dim_source))
    logger.info('The parameter type is {}.'.format(discretization.parameter_type))

    # create training set
    training_set_sampling_strategy = config.get('pymor', 'training_set')
    if training_set_sampling_strategy == 'random':
        num_training_samples = config.getint('pymor', 'num_training_samples')
        training_samples = list(discretization.parameter_space.sample_randomly(num_training_samples))
    else:
        raise ConfigError('unknown \'training_set\' sampling strategy given: \'{}\''.format(training_set_sampling_strategy))

    logger.info('running standard rb for global discretization:')
    reduction_report, data = perform_lrbms(config, discretization, training_samples)
    logger.info(reduction_report)

    # test quality
    logger.info('testing reduction quality:')
    test_set_sampling_strategy = config.get('pymor', 'test_set')
    if test_set_sampling_strategy == 'training':
        test_samples = training_samples
    elif test_set_sampling_strategy == 'random':
        num_test_samples = config.getint('pymor', 'num_test_samples')
        test_samples = list(discretization.parameter_space.sample_randomly(num_test_samples))
    else:
        raise ConfigError('unknown \'test_set\' sampling strategy given: \'{}\''.format(test_set_sampling_strategy))
    test_report = test_quality(config, test_samples, discretization, data, test_set_sampling_strategy)
    logger.info(test_report)
