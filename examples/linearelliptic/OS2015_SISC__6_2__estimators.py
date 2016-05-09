#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

import numpy as np

from pymor.core                 import ImmutableInterface
from pymor.reductors            import reduce_generic_rb
from pymor.playground.reductors import GenericBlockRBReconstructor


class DetailedEstimator(ImmutableInterface):


    def __init__(self, example, wrapper, mu_hat, mu_bar):
        self._example = example
        self._wrapper = wrapper
        self._mu_hat = mu_hat
        self._mu_bar = mu_bar

    def estimate(self, U_h, mu):
        return self._example.estimate(U_h, 'eta_OS2014_*', self._mu_hat, self._mu_bar, mu)


class ReducedEstimator(object):

    def __init__(self, discretization, example, wrapper, mu_hat, mu_bar, norm, compute_ids, return_id):
        self._discretization = discretization
        self._wrapper = wrapper
        self._estimator = DetailedEstimator(example, wrapper, mu_hat, mu_bar)
        self._norm = norm
        self._compute_ids = compute_ids if isinstance(compute_ids, tuple) else (compute_ids,)
        self._return_id = return_id
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
        if 'discretization_error' in self._compute_ids:
            self.add_to_data('discretization_error', mu, example.compute_error(U_h_dune, 'elliptic', mu_dune, mu_bar_dune))
        if 'full_error' in self._compute_ids:
            self.add_to_data('full_error', mu, example.compute_error(U_red_dune, 'elliptic', mu_dune, mu_bar_dune))
        if 'model_reduction_error' in self._compute_ids:
            self.add_to_data('model_reduction_error', mu, self._norm(U_red - U_h)[0])
        # compute estimates
        alpha_mu_mu_bar = example.alpha(mu_dune, mu_bar_dune)
        gamma_mu_mu_bar = example.gamma(mu_dune, mu_bar_dune)
        alpha_mu_mu_hat = example.alpha(mu_dune, mu_hat_dune)
        self.add_to_data('alpha_mu_mu_bar', mu, alpha_mu_mu_bar)
        self.add_to_data('gamma_mu_mu_bar', mu, gamma_mu_mu_bar)
        self.add_to_data('alpha_mu_mu_hat', mu, alpha_mu_mu_hat)
        if 'eta_red' or 'eta_nc_red' in self._compute_ids:
            eta_nc_red = example.estimate(U_red_dune, 'eta_NC_OS2014', mu_hat_dune, mu_bar_dune, mu_dune)
            self.add_to_data('eta_nc_red', mu, eta_nc_red)
        if 'eta_red' or 'eta_r_red' in self._compute_ids:
            eta_r_red  = example.estimate(U_red_dune, 'eta_R_OS2014_*', mu_hat_dune, mu_bar_dune, mu_dune)
            self.add_to_data('eta_r_red',  mu, eta_r_red)
        if 'eta_red' or 'eta_df_red' in self._compute_ids:
            eta_df_red = example.estimate(U_red_dune, 'eta_DF_OS2014_*', mu_hat_dune, mu_bar_dune, mu_dune)
            self.add_to_data('eta_df_red', mu, eta_df_red)
        if 'eta_red' in self._compute_ids:
            eta_red = (1.0/np.sqrt(alpha_mu_mu_bar))*(np.sqrt(gamma_mu_mu_bar)*eta_nc_red
                                                      + eta_r_red
                                                      + (1.0/np.sqrt(alpha_mu_mu_hat))*eta_df_red)
            self.add_to_data('eta_red',    mu, eta_red)
        assert self._return_id in self._compute_ids
        return self.data[self.extension_step][self._return_id][-1]


def reduce_with_estimator(discretization,
                          RB,
                          operator_product=None,
                          vector_product=None,
                          disable_caching=True,
                          extends=None,
                          reduced_estimator=None):
    assert operator_product is None
    rd, _, reduction_data = reduce_generic_rb(discretization,
                                              RB,
                                              vector_product,
                                              disable_caching,
                                              extends)
    rc = GenericBlockRBReconstructor(RB)
    reduced_estimator.extension_step += 1
    reduced_estimator.rc = rc
    rd = rd.with_(estimator=reduced_estimator)
    return rd, rc, reduction_data

