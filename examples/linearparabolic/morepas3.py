#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

import numpy as np

from tempfile import NamedTemporaryFile

from pymor.core.defaults import set_defaults
set_defaults({'pymor.core.cache.default_regions.disk_max_size': 107374182400})

import pymor.core.logger
from pymor.core.logger import getLogger

from simdb.run import new_dataset, add_values, add_logfile

from morepas3__estimate import DetailedAgainstReference, DetailedAgainstWeak
from morepas3__prepare import discretize, prepare

logfile = NamedTemporaryFile(delete=False).name
pymor.core.logger.FILENAME = logfile

logger = getLogger('.morepas3.main')
logger.setLevel('INFO')


config= {'dune_subdomains': '[1 1]',
         'dune_num_elements': '[8 8]',
         'dune_num_refinements': '0',
         'integration_order_time': 2,
         'end_time' : 0.01,
         'nt' : 10,
         'mu_min': 0.1,
         'mu_max': 1,
         'mu_bar': 1,
         'mu_hat': 1,
         'poincare': 1./(np.pi**2),
         'extension_product': 'h1',
         'num_training_samples' : 3,
         'max_rb_size' : 500,
         'target_error': 1e-1,
         'initial_data': 'exp(-((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5))/0.01)'}
new_dataset('morepas3', **config)

logger.info('creating DUNE discretization ...')
detailed_data = prepare(config)
(example, wrapper, initial_data, elliptic_disc, parabolic_disc, bochner_norms,
 mu_min, mu_max, mu_hat, mu_bar) = (
        detailed_data['example'],
        detailed_data['wrapper'],
        detailed_data['initial_data'],
        detailed_data['elliptic_LRBMS_disc'],
        detailed_data['parabolic_disc'],
        detailed_data['bochner_norms'],
        detailed_data['mu_min'],
        detailed_data['mu_max'],
        detailed_data['mu_hat'],
        detailed_data['mu_bar'])
logger.info('  grid has {} subdomains'.format(elliptic_disc.num_subdomains))
logger.info('  discretization has {} DoFs'.format(parabolic_disc.solution_space.dim))
logger.info('  parameter type is {}'.format(parabolic_disc.parameter_type))

# logger.info('visualizing data functions and sample solution ...')
# example.visualize('example')
# U = parabolic_disc.solve(1)
# parabolic_disc.visualize(U, filename='sample_solution')

logger.info('creating reference discretization ...')
config_ref = dict(config)
config_ref['dune_num_elements'] = '[16 16]'
config_ref['nt'] = 2*config['nt']
reference_data = prepare(config_ref)
parabolic_disc_ref, prolongator, bochner_norms_ref = (
        reference_data['parabolic_disc'],
        reference_data['prolongator'],
        reference_data['bochner_norms'])
logger.info('  discretization has {} DoFs'.format(parabolic_disc_ref.solution_space.dim))

logger.info('computing discretization errors ...')
reference_error_computer = DetailedAgainstReference(parabolic_disc_ref,
                                                    prolongator,
                                                    elliptic_disc,
                                                    bochner_norms_ref['elliptic'])

training_samples = list(parabolic_disc.parameter_space.sample_randomly(config['num_training_samples']))
discretization_errors = [reference_error_computer.estimate(parabolic_disc.solve(mu), mu, parabolic_disc)
                         for mu in training_samples]

logger.info('estimating discretization errors ...')
detailed_estimator = DetailedAgainstWeak(example, wrapper, bochner_norms['elliptic'], config['end_time'],
                                         config['poincare'], mu_min, mu_max, mu_hat, mu_max)
discretization_estimates = [detailed_estimator.estimate(parabolic_disc.solve(mu), mu, parabolic_disc)
                            for mu in training_samples]

for ii in np.arange(len(training_samples)):
    print('{}: {} vs. {} ({})'.format(training_samples[ii],
                                      discretization_errors[ii],
                                      discretization_estimates[ii],
                                      discretization_estimates[ii]/discretization_errors[ii]))

# logger.info(' ')

# add_values(time=greedy_data['time'],
#            max_err_mus=greedy_data['max_err_mus'],
#            extensions=greedy_data['extensions'],
#            max_errs=greedy_data['max_errs'],
#            basis_sizes=[len(local_RB) for local_RB in RB])
# add_logfile(logfile)

