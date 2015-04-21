#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

import numpy as np

import pymor.core
from pymor.core.logger import getLogger
from pymor.operators.constructions import induced_norm
from pymor.operators.constructions import LincombOperator
from pymor.parameters.spaces import CubicParameterSpace

from dune.pymor.core import wrap_module

import linearellipticexampleOS2014 as dune_module


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

    example.visualize(cfg['dune_example'])

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

