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
from pymor.core import getLogger, ImmutableInterface
from pymor.la import induced_norm
from pymor.parameters import CubicParameterSpace, Parametric
from pymor.operators import LincombOperator
from pymor.reductors import reduce_generic_rb
from pymor.playground.la import BlockVectorArray
from pymor.playground.algorithms import gram_schmidt_block_basis_extension, trivial_block_basis_extension
from pymor.playground.reductors import GenericBlockRBReconstructor

from dune.pymor.core import wrap_module

from simdb.run import new_dataset, add_values, add_data, add_logfile

import linearellipticexampleOS2014 as dune_module


dune_config = {'dune_partitioning': '[1 1 1]',
               'dune_num_refinements': 2,
               'dune_products': ['l2', 'h1_semi', 'elliptic'],
               'dune_log_info_level': -1,
               'dune_log_debug_level': -1,
               'dune_log_enable_warnings': True,
               'dune_linear_solver_options': {'type': 'bicgstab.ilut', 'precision': '1e-14'},
               'dune_example': 'OS2014Example'}
config = {'num_training_samples': 10,
          'mu_hat_value': 0.1,
          'mu_bar_value': 0.1,
          'greedy_use_estimator': True,
          'greedy_target_error': 1e-14,
          'initialize_with_one': True,
          'estimator_return': 'eta_red'}
DATASET_ID = dune_config['dune_example'] + '_oversampling_test'


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


class ReducedEstimator(object):

    def __init__(self, discretization, example, wrapper, mu_hat, mu_bar, norm):
        self._discretization = discretization
        self._wrapper = wrapper
        self._estimator = DetailedEstimator(example, wrapper, mu_hat, mu_bar)
        self._norm = norm
        self.extension_step = -1
        self.rc = None
        self.data = {}

    def add_to_data(self, key, mu, value):
        if not self.data.has_key(self.extension_step):
            self.data[self.extension_step] = {}
        if not self.data[self.extension_step].has_key(key):
            self.data[self.extension_step][key] = []
        if not self.data[self.extension_step].has_key(key + '_mus'):
            self.data[self.extension_step][key + '_mus'] = []
        self.data[self.extension_step][key].append(value)
        self.data[self.extension_step][key + '_mus'].append(mu)

    def estimate(self, U, mu, discretization):
        U_red = self.rc.reconstruct(U)
        assert len(U_red) == 1
        U_red_global = self._discretization.globalize_vectors(U_red)
        U_red_dune = U_red_global._list[0]._impl
        U_h = self._discretization.solve(mu)
        U_h_global = self._discretization.globalize_vectors(U_h)
        assert len(U_h_global) == 1
        U_h_dune = U_h_global._list[0]._impl
        # compute errors
        example = self._estimator._example
        mu_dune = self._wrapper.dune_parameter(mu)
        mu_bar_dune = self._estimator._mu_bar
        mu_hat_dune = self._estimator._mu_hat
        self.add_to_data('discretization_error', mu, example.compute_error(U_h_dune, 'elliptic', mu_dune, mu_bar_dune))
        self.add_to_data('full_error', mu, example.compute_error(U_red_dune, 'elliptic', mu_dune, mu_bar_dune))
        self.add_to_data('model_reduction_error', mu, self._norm(U_red - U_h)[0])
        # compute estimates
        alpha_mu_mu_bar = example.alpha(mu_dune, mu_bar_dune)
        gamma_mu_mu_bar = example.gamma(mu_dune, mu_bar_dune)
        alpha_mu_mu_hat = example.alpha(mu_dune, mu_hat_dune)
        self.add_to_data('alpha_mu_mu_bar', mu, alpha_mu_mu_bar)
        self.add_to_data('gamma_mu_mu_bar', mu, gamma_mu_mu_bar)
        self.add_to_data('alpha_mu_mu_hat', mu, alpha_mu_mu_hat)
        eta_nc_red = example.estimate(U_red_dune, 'eta_NC_OS2014', mu_hat_dune, mu_bar_dune, mu_dune)
        eta_r_red  = example.estimate(U_red_dune, 'eta_R_OS2014_*', mu_hat_dune, mu_bar_dune, mu_dune)
        eta_df_red = example.estimate(U_red_dune, 'eta_DF_OS2014_*', mu_hat_dune, mu_bar_dune, mu_dune)
        eta_red = (1.0/np.sqrt(alpha_mu_mu_bar))*(np.sqrt(gamma_mu_mu_bar)*eta_nc_red
                                                  + eta_r_red
                                                  + (1.0/np.sqrt(alpha_mu_mu_hat))*eta_df_red)
        self.add_to_data('eta_nc_red', mu, eta_nc_red)
        self.add_to_data('eta_r_red',  mu, eta_r_red)
        self.add_to_data('eta_df_red', mu, eta_df_red)
        self.add_to_data('eta_red',    mu, eta_red)
        # self.add_to_data('eta_red', mu, example.estimate(U_red_dune, 'eta_OS2014_*', mu_hat_dune, mu_bar_dune, mu_dune))
        return self.data[self.extension_step][config['estimator_return']][-1]


def reduce_with_estimator(discretization,
                          RB,
                          operator_product=None,
                          vector_product=None,
                          disable_caching=True,
                          extends=None,
                          reduced_estimator=None):
    rd, _, reduction_data = reduce_generic_rb(discretization,
                                               RB,
                                               operator_product,
                                               vector_product,
                                               disable_caching,
                                               extends)
    rc = GenericBlockRBReconstructor(RB)
    reduced_estimator.extension_step += 1
    reduced_estimator.rc = rc
    rd = rd.with_(estimator=reduced_estimator)
    return rd, rc, reduction_data


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
                         max_extensions=len(training_samples)+1,
                         target_error=cfg['greedy_target_error'])
    add_values(time=greedy_data['time'],
               max_err=greedy_data['max_err'],
               max_err_mus=greedy_data['max_err_mus'],
               extensions=greedy_data['extensions'],
               max_err_mu=greedy_data['max_err_mu'],
               max_errs=greedy_data['max_errs'],
               estimator_data=reduced_estimator.data)

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

