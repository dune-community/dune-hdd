#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from tempfile import NamedTemporaryFile

import pymor.core
from pymor.core import getLogger

from simdb.run import new_dataset, add_logfile

from OS2014_prepare import prepare
from OS2014_offline import offline_phase
from OS2014_online  import online_phase


config = {'dune_partitioning': '[25 5 1]',
          'dune_num_refinements': 0,
          'dune_oversampling_layers': 8,
          'dune_products': ['elliptic', 'penalty'],
          'dune_log_info_level': -1,
          'dune_log_debug_level': -1,
          'dune_log_enable_warnings': True,
          'dune_linear_solver_options': {'type': 'bicgstab.ilut', 'precision': '1e-14'},
          'dune_example': 'Spe10Model1Example',
          'parameter_range': (0.1, 1.0),
          'compute_some_solution_norms': False,
          'initialize_basis_with': 1,
          'training_sampling_strategy': 'uniform',
          'num_training_samples': 2,
          'greedy_max_extensions': 2,
          'extension_product': ('elliptic', 'penalty'),
          'greedy_error_norm': 'elliptic',
          'mu_hat_value': 0.1,
          'mu_bar_value': 0.1,
          'greedy_target_error': 1e-10,
          'greedy_use_estimator': True,
          'estimator_compute': 'eta_red',
          'estimator_return': 'eta_red',
          'num_test_samples': 10,
          'estimate_some_errors': True,
          'local_indicators': 'eta_red',
          'marking_strategy': 'doerfler_and_age',
          'marking_max_age': 3,
          'doerfler_marking_theta': 0.5,
          'local_boundary_values': 'dirichlet',
          'online_target_error': 5,
          'online_max_extensions': 50}
DATASET_ID = config['dune_example'] + '_online_enrichment_test'

pymor.core.logger.MAX_HIERACHY_LEVEL = 2
getLogger('pymor.WrappedDiscretization').setLevel('WARN')
getLogger('pymor.algorithms').setLevel('INFO')
getLogger('dune.pymor.discretizations').setLevel('WARN')


if __name__ == '__main__':

    logfile = NamedTemporaryFile(delete=False).name
    pymor.core.logger.FILENAME = logfile
    new_dataset(DATASET_ID, **config)

    detailed_data = prepare(config)
    print('')
    offline_data  = offline_phase(config, detailed_data)
    print('')
    _             = online_phase(config, detailed_data, offline_data)

    # this should be the last action (to really capture all logs)
    add_logfile(logfile)

