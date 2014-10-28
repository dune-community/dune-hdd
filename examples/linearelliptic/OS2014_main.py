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

from pymor.algorithms import greedy, gram_schmidt_basis_extension
import pymor.core
from pymor.core import getLogger, ImmutableInterface
from pymor.la import induced_norm
from pymor.parameters import CubicParameterSpace, Parametric
from pymor.operators import LincombOperator
from pymor.reductors import reduce_generic_rb

from dune.pymor.core import wrap_module

from simdb.run import new_dataset, add_values, add_data, add_logfile

import linearellipticexampleOS2014 as dune_module


dune_config = {'dune_partitioning': '[1 1 1]',
               'dune_num_refinements': 2,
               'dune_products': ['l2', 'h1_semi', 'elliptic', 'energy'],
               'dune_log_info_level': -1,
               'dune_log_debug_level': -1,
               'dune_log_enable_warnings': True}
config = {'num_training_samples': 100,
          'mu_bar_value': 1.0,
          'greedy_use_estimator': False,
          'greedy_target_error': 1e-12,
          'initialize_with_one': False}


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
    example = dune_module.OS2014Example(partitioning=cfg['dune_partitioning'],
                                        num_refinements=cfg['dune_num_refinements'],
                                        products=cfg['dune_products'],
                                        info_log_levels=cfg['dune_log_info_level'],
                                        debug_log_levels=cfg['dune_log_debug_level'],
                                        enable_warnings=cfg['dune_log_enable_warnings'])
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




def run_experiment(example, wrapper, cfg, product, norm):

    logger, logfile = get_logger()


    discretization = wrapper[example.discretization()]
    discretization = discretization.with_(
                    parameter_space=CubicParameterSpace(discretization.parameter_type, 0.1, 1.0))
    logger.info('the discretization has {} DoFs.'.format(discretization.solution_space.dim))
    logger.info('the parameter type is {}.'.format(discretization.parameter_type))
    logger.info('')

    logger.info('computing solution norms:')
    mu_bar_dune = wrapper.DuneParameter('mu', cfg['mu_bar_value'])
    products = {}
    for nm, pr in discretization.products.items():
        if not pr.parametric:
            products[nm] = pr
        else:
            products[nm] = wrapper[pr._impl.freeze_parameter(mu_bar_dune)]

    product, product_name = create_product(products, product)
    norm, norm_name = create_product(products, norm)
    norm = induced_norm(norm)
    cfg['extension_product'] = product_name
    cfg['greedy_error_norm'] = norm_name

    new_dataset('greedy', 'OS2014', **cfg)

    training_samples = list(discretization.parameter_space.sample_randomly(cfg['num_training_samples']))
    solution_norms = [norm(discretization.solve(mu)) for mu in training_samples]
    logger.info('  range:              [{}, {}]'.format(np.amin(solution_norms), np.amax(solution_norms)))
    logger.info('  median:              {}'.format(np.median(solution_norms)))
    logger.info('  mean:                {}'.format(np.mean(solution_norms)))
    logger.info('  standard deviation:  {}'.format(np.std(solution_norms)))
    add_values(training_samples=training_samples,
               solution_norms=solution_norms)
    logger.info('')

    extension_algorithm=partial(gram_schmidt_basis_extension, product=product)
    initial_basis = discretization.functionals['rhs'].source.empty()
    if cfg['initialize_with_one']:
        one = wrapper.vector_array(wrapper[example.project('1')])
        initial_basis = extension_algorithm(initial_basis, one)[0]

    greedy_data = greedy(discretization,
                         reduce_generic_rb,
                         training_samples,
                         initial_basis=initial_basis,
                         use_estimator=cfg['greedy_use_estimator'],
                         error_norm=norm,
                         max_extensions=len(training_samples)+1,
                         target_error=cfg['greedy_target_error'])
    add_values(time=greedy_data['time'],
               max_err=greedy_data['max_err'],
               max_err_mus=greedy_data['max_err_mus'],
               extensions=greedy_data['extensions'],
               max_err_mu=greedy_data['max_err_mu'],
               max_errs=greedy_data['max_errs'])

    # this should be the last action
    add_logfile(logfile)


if __name__ == '__main__':

    print('initializing dune module... ')
    example, wrapper = init_dune(dune_config)

    for product, norm in {'l2': 'l2',
                          ('l2', 'h1_semi'): ('l2', 'h1_semi'),
                          'energy': 'energy'}.items():
        cfg = config.copy()
        for kk, vv in dune_config.items():
            assert not cfg.has_key(kk)
            cfg[kk] = vv

        print('')
        print('running experiment with product \'{}\' and norm \'{}\':'.format(product, norm))
        run_experiment(example, wrapper, cfg, product, norm)

