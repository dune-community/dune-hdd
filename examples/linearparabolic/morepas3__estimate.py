#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from functools import partial

import numpy as np

from pymor.grids.oned import OnedGrid

from dune.pymor.la.container import make_listvectorarray


def L2_time_L2_space_partial_t_estimate(l2_product, T, U):
    result = 0.
    time_grid = OnedGrid(domain=(0., T), num_intervals=len(U)-1)
    dt = time_grid.volumes(0)[0]
    for n in np.arange(1, time_grid.size(1)):
        diff = make_listvectorarray(U._list[n] - U._list[n - 1])
        result += l2_product.apply2(diff, diff)[0][0]
    result /= dt
    return np.sqrt(result)


def time_residual_estimate(wrapper, T, b, l2_product, U, mu):
    result = 0.
    time_grid = OnedGrid(domain=(0., T), num_intervals=len(U)-1)
    dt = time_grid.volumes(0)[0]
    for n in np.arange(1, time_grid.size(1)):
        b_func = b.apply(make_listvectorarray(U._list[n] - U._list[n - 1]), mu=mu)
        b_riesz = l2_product.apply_inverse(b_func)
        result += l2_product.apply2(b_riesz, b_riesz)[0][0]
    result *= dt/3.
    return np.sqrt(result)


def elliptic_reconstruction_estimate(example, T, U, mu_min, mu_max, mu_tilde, mu):
    result = 0.
    time_grid = OnedGrid(domain=(0., T), num_intervals=len(U)-1)
    dt = time_grid.volumes(0)[0]
    for n in np.arange(time_grid.size(1)):
        result += example.elliptic_reconstruction_estimate(U._list[n]._impl, mu, mu, mu_tilde, mu, mu, '')**2
        # result += example.elliptic_reconstruction_estimate(U._list[n]._impl, mu_min, mu_max, mu_tilde, mu, mu)**2
    result *= dt/3.
    return 2*np.sqrt(result)


class DetailedAgainstReference(object):

    def __init__(self, reference_disc, prolong, coarse_elliptic_disc, bochner_norm):
        self._reference_disc = reference_disc
        self._prolong = prolong
        self._coarse_elliptic_disc = coarse_elliptic_disc
        self._bochner_norm = bochner_norm

    def estimate(self, U, mu, disc):
        U_reference = self._reference_disc.solve(mu)
        U_coarse_prologated = self._prolong(self._coarse_elliptic_disc, U)
        return self._bochner_norm(U_reference - U_coarse_prologated)


class DetailedAgainstWeak(object):

    def __init__(self, example, wrapper, L2_elliptic_bochner_norm, l2_product, T, mu_min, mu_max, mu_hat, mu_bar, mu_tilde):
        self._example = example
        self._wrapper = wrapper
        self._L2_elliptic_bochner_norm = L2_elliptic_bochner_norm
        self._l2_product = l2_product
        self._T = T
        self._mu_min = wrapper.dune_parameter(mu_min)
        self._mu_max = wrapper.dune_parameter(mu_max)
        self._mu_hat = wrapper.dune_parameter(mu_hat)
        self._mu_bar = wrapper.dune_parameter(mu_bar)
        self._mu_tilde = wrapper.dune_parameter(mu_tilde)

    def compute_indicators(self, U, mu, disc):
        p_N = U
        _mu = mu
        mu = self._wrapper.dune_parameter(mu)
        p_N_c = make_listvectorarray([self._wrapper[self._example.oswald_interpolate(p._impl)] for p in p_N._list])
        p_N_d = p_N - p_N_c

        c_eps_mu_hat = self._example.min_diffusion_ev(self._mu_hat)
        alpha_mu_mu_hat = self._example.alpha(mu, self._mu_hat)
        C_P_Omega = 2*self._example.domain_diameter()

        e_c_0_norm = 0.
        dt_p_N_d_norm = L2_time_L2_space_partial_t_estimate(self._l2_product, self._T, p_N_d)
        eps_norm = elliptic_reconstruction_estimate(self._example,
                                                    self._T,
                                                    p_N,
                                                    self._mu_min, self._mu_max, self._mu_tilde, mu)
        p_N_d_norm = self._L2_elliptic_bochner_norm(p_N_d, mu=_mu)
        R_T_norm = time_residual_estimate(self._wrapper,
                                          self._T,
                                          disc.operator,
                                          disc.l2_product,
                                          p_N,
                                          _mu)

        return {'c_eps_mu_hat': c_eps_mu_hat,
                'alpha_mu_mu_hat': alpha_mu_mu_hat,
                'C_P_Omega': C_P_Omega,
                'e_c_0_norm': e_c_0_norm,
                'dt_p_N_d_norm': dt_p_N_d_norm,
                'eps_norm': eps_norm,
                'p_N_d_norm': p_N_d_norm,
                'R_T_norm': R_T_norm}

    def estimate(self, U, mu, disc):
        alpha_mu_mu_bar = self._example.alpha(self._wrapper.dune_parameter(mu), self._mu_bar)
        indicators = self.compute_indicators(U, mu, disc)
        (e_c_0_norm, C_P_Omega, c_eps_mu_hat, alpha_mu_mu_hat, dt_p_N_d_norm, eps_norm, p_N_d_norm, R_T_norm) = (
                indicators['e_c_0_norm'], indicators['C_P_Omega'], indicators['c_eps_mu_hat'],
                indicators['alpha_mu_mu_hat'], indicators['dt_p_N_d_norm'], indicators['eps_norm'],
                indicators['p_N_d_norm'], indicators['R_T_norm'])
        return  1./np.sqrt(alpha_mu_mu_bar) * (  e_c_0_norm
                                               + (2.*C_P_Omega*c_eps_mu_hat)/alpha_mu_mu_hat * dt_p_N_d_norm
                                               + (np.sqrt(5.) + 1.) * eps_norm
                                               + np.sqrt(5.) * p_N_d_norm
                                               + 2./alpha_mu_mu_hat * R_T_norm)

