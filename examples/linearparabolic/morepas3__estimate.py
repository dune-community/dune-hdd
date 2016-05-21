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


def dual_energy_space_norm_squared_estimate(l2_product, C_P, alpha, c_eps, U):
    U = make_listvectorarray(U)
    return (C_P/(alpha * c_eps))**2 * l2_product.apply2(U, U)[0][0]


def L2_time_dual_energy_space_partial_t_estimate(dual_energy2_estimate, T, U):
    result = 0.
    time_grid = OnedGrid(domain=(0., T), num_intervals=len(U)-1)
    for n in np.arange(1, time_grid.size(1)):
        dt = time_grid.volumes(0)[n - 1]
        result += 1./dt * dual_energy2_estimate(U._list[n] - U._list[n - 1])
    return np.sqrt(result)


def time_residual_estimate(wrapper, T, B, l2_product, U, mu):
    B =  wrapper[B._impl.freeze_parameter(mu)]
    result = 0.
    time_grid = OnedGrid(domain=(0., T), num_intervals=len(U)-1)
    for n in np.arange(1, time_grid.size(1)):
        B_func = B.apply(make_listvectorarray(U._list[n] - U._list[n - 1]))
        B_riesz = l2_product.apply_inverse(B_func)
        dt = time_grid.volumes(0)[n - 1]
        result += dt/3. * l2_product.apply2(B_riesz, B_riesz)[0][0]
    return np.sqrt(result)


def elliptic_reconstruction_estimate(example, T, U, mu_min, mu_max, mu_hat, mu_bar, mu):
    result = 0.
    time_grid = OnedGrid(domain=(0., T), num_intervals=len(U)-1)
    dt = time_grid.volumes(0)[0]
    for n in np.arange(time_grid.size(1)):
        result += dt/3. * example.elliptic_reconstruction_estimate(U._list[n]._impl,
                                                                   mu_min, mu_max, mu_hat, mu_bar,
                                                                   mu)**2
    return np.sqrt(result)


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

    def __init__(self, example, wrapper, L2_elliptic_bochner_norm, T, C_P, mu_min, mu_max, mu_hat, mu_bar):
        self._example = example
        self._wrapper = wrapper
        self._L2_elliptic_bochner_norm = L2_elliptic_bochner_norm
        self._T = T
        self._C_P = C_P
        self._mu_min = wrapper.dune_parameter(mu_min)
        self._mu_max = wrapper.dune_parameter(mu_max)
        self._mu_hat = wrapper.dune_parameter(mu_hat)
        self._mu_bar = wrapper.dune_parameter(mu_bar)

    def estimate(self, U, mu, disc):
        p_N = U
        mu = self._wrapper.dune_parameter(mu)
        p_N_c = make_listvectorarray([self._wrapper[self._example.oswald_interpolate(p._impl)] for p in p_N._list])
        p_N_d = p_N - p_N_c

        c_eps_mu_hat = self._example.min_diffusion_ev(self._mu_hat)
        C_eps_mu_hat = self._example.max_diffusion_ev(self._mu_hat)
        alpha_mu_mu_hat = self._example.alpha(mu, self._mu_hat)
        gamma_mu_mu_hat = self._example.gamma(mu, self._mu_hat)

        C_b = 1. #gamma_mu_mu_hat * C_eps_mu_hat
        C = np.sqrt(3*C_b**2 + 2)
        C_V = self._C_P/(alpha_mu_mu_hat * c_eps_mu_hat)

        e_c_0_norm = 0.
        dt_p_N_d_norm = L2_time_dual_energy_space_partial_t_estimate(partial(dual_energy_space_norm_squared_estimate,
                                                                             disc.l2_product,
                                                                             self._C_P,
                                                                             alpha_mu_mu_hat,
                                                                             c_eps_mu_hat),
                                                                     self._T,
                                                                     p_N_d)
        eps_norm = elliptic_reconstruction_estimate(self._example,
                                                    self._T,
                                                    p_N,
                                                    self._mu_min, self._mu_max, self._mu_hat, self._mu_bar, mu)
        p_N_d_norm = self._L2_elliptic_bochner_norm(p_N_d)
        R_T_norm = time_residual_estimate(self._wrapper,
                                          self._T,
                                          disc.operator,
                                          disc.l2_product,
                                          p_N,
                                          mu)

        return 2.*dt_p_N_d_norm + (C + 1)*eps_norm + C*p_N_d_norm + 2.*C_V*R_T_norm

