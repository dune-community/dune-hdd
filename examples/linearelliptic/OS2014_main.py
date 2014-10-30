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
               'dune_products': ['l2', 'h1_semi', 'elliptic', 'energy', 'penalty'],
               'dune_log_info_level': -1,
               'dune_log_debug_level': -1,
               'dune_log_enable_warnings': True,
               'dune_linear_solver_options': {'type': 'bicgstab.ilut', 'precision': '1e-14'}}
config = {'num_training_samples': 100,
          'mu_hat_value': 1.0,
          'mu_bar_value': 1.0,
          'greedy_use_estimator': False,
          'greedy_target_error': 1e-14,
          'initialize_with_one': True}


pymor.core.logger.MAX_HIERACHY_LEVEL = 2
getLogger('pymor.WrappedDiscretization').setLevel('WARN')
getLogger('pymor.algorithms').setLevel('INFO')
getLogger('dune.pymor.discretizations').setLevel('WARN')


class DetailedEstimator(ImmutableInterface):


    def __init__(self, example, wrapper, mu_hat, mu_bar):
        self._example = example
        self._wrapper = wrapper
        self._mu_hat = mu_hat
        self._mu_bar = mu_bar

    def estimate(self, U_h, mu):
        return self._example.estimate(U_h, 'eta_OS2014_*', self._mu_hat, self._mu_bar, mu)


class ReducedEstimator(ImmutableInterface):

    def __init__(self, discretization, example, wrapper, mu_hat, mu_bar, norm, rc):
        self._discretization = discretization
        self._wrapper = wrapper
        self._estimator = DetailedEstimator(example, wrapper, mu_hat, mu_bar)
        self._norm = norm
        self._rc = rc

    def estimate(self, U, mu, discretization):
        mu_dune = self._wrapper.dune_parameter(mu)
        # mu_bar_dune = self._estimator._mu_bar
        # mu_hat_dune = self._estimator._mu_hat
        # example = self._estimator._example
        # print('mu     = {}'.format(mu_dune.report()))
        # print('mu_bar = {}'.format(mu_bar_dune.report()))
        # print('mu_hat = {}'.format(mu_hat_dune.report()))
        U_red = self._rc.reconstruct(U)
        assert len(U_red) == 1
        U_red_dune = U_red._list[0]._impl
        # U_h = self._discretization.solve(mu)
        # assert len(U_h) == 1
        # U_h_dune = U_h._list[0]._impl
        # print('    || U(mu) - U_red(mu) ||_mu_bar = {}'.format(example.compute_error(U_red_dune,
        #                                                                              'elliptic',
        #                                                                              mu_dune,
        #                                                                              mu_bar_dune)))
        # print('    || U(mu) - U_h(mu)   ||_mu_bar = {}'.format(example.compute_error(U_h_dune,
        #                                                                              'elliptic',
        #                                                                              mu_dune,
        #                                                                              mu_bar_dune)))
        return self._estimator.estimate(U_red_dune, mu_dune)
        # print('         eta[mu_bar,mu_hat](U_red) = {}'.format(eta))
        # if self._norm is not None:
        #     err = self._norm(U_red - U_h)[0]
        #     print('  || U_h(mu) - U_red(mu) ||_mu_bar = {}'.format(err))
        #     return err
        # else:
        #     return eta


def reduce_with_estimator(discretization,
                          RB,
                          operator_product=None,
                          vector_product=None,
                          disable_caching=True,
                          extends=None,
                          detailed_discretization=None,
                          example=None,
                          wrapper=None,
                          mu_hat=None,
                          mu_bar=None,
                          norm=None):
    rd, rc, reduction_data = reduce_generic_rb(discretization,
                                               RB,
                                               operator_product,
                                               vector_product,
                                               disable_caching,
                                               extends)
    rd = rd.with_(estimator=ReducedEstimator(detailed_discretization, example, wrapper, mu_hat, mu_bar, norm, rc))
    return rd, rc, reduction_data


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

    discretization = example.discretization()
    linear_solver_options = cfg['dune_linear_solver_options']
    dune_linear_solver_options = discretization.solver_options(linear_solver_options['type'])
    for kk, vv in linear_solver_options.items():
        dune_linear_solver_options.set(kk, vv, True)
    return wrapper[discretization].with_(solver_options=dune_linear_solver_options)


def run_experiment(example, wrapper, cfg, product, norm):

    logger, logfile = get_logger()

    discretization = create_discretization(example, wrapper, cfg)
    discretization = discretization.with_(
                    parameter_space=CubicParameterSpace(discretization.parameter_type, 0.1, 1.0))
    logger.info('the discretization has {} DoFs.'.format(discretization.solution_space.dim))
    logger.info('the parameter type is {}.'.format(discretization.parameter_type))
    logger.info('')

    logger.info('computing solution norms:')
    mu_hat_dune = wrapper.DuneParameter('mu', cfg['mu_hat_value'])
    mu_bar_dune = wrapper.DuneParameter('mu', cfg['mu_bar_value'])
    products = {}
    for nm, pr in discretization.products.items():
        if not pr.parametric:
            products[nm] = pr
        else:
            products[nm] = wrapper[pr._impl.freeze_parameter(mu_bar_dune)]

    product, product_name = create_product(products, product)
    norm, norm_name = create_product(products, norm)
    norm = induced_norm(norm) if norm is not None else None
    cfg['extension_product'] = product_name
    cfg['greedy_error_norm'] = norm_name

    new_dataset('greedy', 'OS2014', **cfg)

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

    extension_algorithm=partial(gram_schmidt_basis_extension, product=product)
    initial_basis = discretization.functionals['rhs'].source.empty()
    if cfg['initialize_with_one']:
        one = wrapper.vector_array(wrapper[example.project('1')])
        initial_basis = extension_algorithm(initial_basis, one)[0]

    greedy_data = greedy(discretization,
                         partial(reduce_with_estimator,
                                 detailed_discretization=discretization,
                                 example=example,
                                 wrapper=wrapper,
                                 mu_hat=mu_hat_dune,
                                 mu_bar=mu_bar_dune,
                                 norm=norm),
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

    # this should be the last action (to really capture all logs)
    add_logfile(logfile)


if __name__ == '__main__':

    print('initializing dune module... ')
    example, wrapper = init_dune(dune_config)

    for product, norm in ((None, None),
                          ('l2', 'l2'),
                          (('l2', 'h1_semi'), ('l2', 'h1_semi')),
                          ('energy', 'energy'),
                          (('elliptic', 'penalty'), 'elliptic')):
        cfg = config.copy()
        for kk, vv in dune_config.items():
            assert not cfg.has_key(kk)
            cfg[kk] = vv

        print('')
        print('running experiment with product \'{}\' and norm \'{}\':'.format(product, norm))
        run_experiment(example, wrapper, cfg, product, norm)
