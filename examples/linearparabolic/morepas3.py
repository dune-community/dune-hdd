#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

import numpy as np

from functools import partial
from tempfile import NamedTemporaryFile

from simdb.run import new_dataset, add_values

from pymor.core.defaults import set_defaults
set_defaults({'pymor.core.cache.default_regions.disk_max_size': 107374182400})

import pymor.core.logger
from pymor.algorithms.error import reduction_error_analysis
from pymor.playground.algorithms.blockbasisextension import gram_schmidt_block_basis_extension
from pymor.core.logger import getLogger
from pymor.vectorarrays.block import BlockVectorArray

from dune.pymor.la.container import make_listvectorarray

from morepas3__estimate import DetailedAgainstReference, DetailedAgainstWeak
from morepas3__prepare import discretize, prepare
from morepas3__reduce import reduce_pod_greedy

logfile = NamedTemporaryFile(delete=False).name
pymor.core.logger.FILENAME = logfile

for logger_id in ('.morepas3.main', '.morepas3.estimate'):
    getLogger(logger_id).setLevel('INFO')
for logger_id in ('pymor.algorithms.gram_schmidt', 'pymor.algorithms.pod'):
    getLogger(logger_id).setLevel('WARN')
logger = getLogger('.morepas3.main')

config= {'dune_num_elements': '[100 20]',
         'dune_num_partitions': '[5 1]',
         'integration_order_time': 2,
         'end_time' : 0.05,
         'nt' : 10,
         'mu_min': 0.1,
         'mu_max': 1.,
         'mu_bar': 0.1,
         'mu_hat': 0.1,
         'mu_tilde': 0.1,
         'extension_product': 'h1',
         'initial_basis': 0,
         'num_training_samples' : 10,
         'num_test_samples' : 10,
         'max_rb_size' : 999,
         'target_error': 1,
         'initial_data': '0'}
new_dataset('morepas3', **config)


# def eoc_study():
#     refinements = (('[4 4]', 10),
#                    ('[8 8]', 20),
#                    ('[16 16]', 40),
#                    ('[32 32]', 80),
#                    ('[64 64]', 160),
#                    ('[128 128]', 320))

#     def level_disc(init, dx_nt):
#         level_cfg = dict(config)
#         level_cfg['dune_num_elements'] = dx_nt[0]
#         level_cfg['nt'] = dx_nt[1]
#         level_cfg['initial_data'] = init
#         return prepare(level_cfg)

#     def compute_error(level_data, reference_data, mu):
#         reference_error_computer = DetailedAgainstReference(reference_data['parabolic_disc'],
#                                                             reference_data['prolongator'],
#                                                             level_data['elliptic_disc'],
#                                                             partial(reference_data['bochner_norms']['elliptic_penalty'], mu=mu))
#         return reference_error_computer.estimate(level_data['parabolic_disc'].solve(mu), mu, level_data['parabolic_disc'])

#     def estimate_error(level_data, mu):
#         detailed_estimator = DetailedAgainstWeak(level_data['example'], level_data['wrapper'],
#                                                  level_data['bochner_norms']['elliptic_penalty'], level_data['space_products']['l2'],
#                                                  config['end_time'],
#                                                  level_data['mu_min'], level_data['mu_max'], level_data['mu_hat'],
#                                                  level_data['mu_bar'], level_data['mu_tilde'])
#         return detailed_estimator.compute_indicators(level_data['parabolic_disc'].solve(mu), mu, level_data['parabolic_disc'])

#     logger.info('creating DUNE discretization ...')
#     level_data = [0, 1, 2, 3]
#     level_data[0] = level_disc(config['initial_data'], refinements[0])
#     for ii in range(1, len(level_data)):
#         level_data[ii] = level_disc((level_data[0]['elliptic_disc'], level_data[0]['initial_data']), refinements[ii])
#     reference_data = level_disc((level_data[0]['elliptic_disc'], level_data[0]['initial_data']), refinements[-1])
#     logger.info('  reference discretization has {} DoFs'.format(reference_data['parabolic_disc'].solution_space.dim))

#     logger.info('computing errors ...')
#     mu = level_data[0]['parabolic_disc'].parse_parameter(1)
#     errors = [compute_error(lvl_data, reference_data, mu) for lvl_data in level_data]

#     logger.info('computing estimates ...')
#     estimates = [estimate_error(lvl_data, mu) for lvl_data in level_data]

#     quantities = {}
#     for kk in estimates[0].keys():
#         quantities[kk] = [None for ii in np.arange(len(estimates))]
#     quantities['error'] = [None for ii in np.arange(len(estimates))]
#     for lvl in np.arange(len(estimates)):
#         lvl_data = estimates[lvl]
#         for kk, vv in lvl_data.items():
#             quantities[kk][lvl] = vv
#         quantities['error'][lvl] = errors[lvl]


#     logger.info('  finished')

#     grid_widths = [lvl_data['example'].max_grid_width() for lvl_data in level_data]

#     def eoc(values):
#         eocs = [None for ii in np.arange(len(values))]
#         for ii in np.arange(1, len(values)):
#             try:
#                 eocs[ii] = np.log(values[ii]/values[ii - 1]) / np.log(grid_widths[ii]/grid_widths[ii - 1])
#             except ZeroDivisionError:
#                 pass
#         return eocs

#     logger.info('errors and estimates:')
#     for kk, vv in quantities.items():
#         print('{}: {}'.format(kk, vv))

#     logger.info(' ')
#     logger.info('EOCs:')
#     for kk, vv in quantities.items():
#         print('{}: {}'.format(kk, eoc(vv)[1:]))


logger.info('creating DUNE discretization ...')
detailed_data = prepare(config)
(example, wrapper, initial_data, elliptic_disc, parabolic_disc, bochner_norms, space_products,
 mu_min, mu_max, mu_hat, mu_bar, mu_tilde) = (
        detailed_data['example'], detailed_data['wrapper'], detailed_data['initial_data'],
        detailed_data['elliptic_LRBMS_disc'], detailed_data['parabolic_disc'], detailed_data['bochner_norms'],
        detailed_data['space_products'], detailed_data['mu_min'], detailed_data['mu_max'], detailed_data['mu_hat'],
        detailed_data['mu_bar'], detailed_data['mu_tilde'])
logger.info('  grid has {} subdomains'.format(elliptic_disc.num_subdomains))
logger.info('  discretization has {} DoFs'.format(parabolic_disc.solution_space.dim))
logger.info('  parameter type is {}'.format(parabolic_disc.parameter_type))

training_samples = list(parabolic_disc.parameter_space.sample_uniformly(config['num_training_samples']))

# logger.info('visualizing data functions and sample solution ...')
# example.visualize('example')
# for mu in (config['mu_min'], config['mu_max']):
#     mu_str = str(mu)
#     mu = parabolic_disc.parse_parameter(mu)
#     U = parabolic_disc.solve(mu)
#     parabolic_disc.visualize(U, filename='sample_solution_{}'.format(mu_str))

logger.info('estimating discretization errors ...')
detailed_data['estimator'] = DetailedAgainstWeak(example, wrapper, bochner_norms['elliptic_penalty'], space_products['l2'],
                                                 config['end_time'], mu_min, mu_max, mu_hat, mu_bar, mu_tilde)
estimates = [detailed_data['estimator'].estimate(parabolic_disc.solve(parabolic_disc.parse_parameter(mu)),
                                                 parabolic_disc.parse_parameter(mu),
                                                 parabolic_disc)
             for mu in (config['mu_min'], config['mu_max'])]
logger.info('  range: {} to {}'.format(np.min(estimates), np.max(estimates)))
if np.max(estimates) > config['target_error']:
    logger.warn(('target error for greedy ({}) is below max estimated discretization error ({}), '
                 + 'greedy will most likely fail!').format(config['target_error'], np.max(estimates)))

if config['initial_basis'] < 0:
    config['initial_basis'] = 0
elif config['initial_basis'] > 1:
    config['initial_basis'] = 1
logger.info('initializing local reduced bases of order up to {} ...'.format(config['initial_basis']))
local_products = [elliptic_disc.local_product(ss, config['extension_product'])
                  for ss in np.arange(elliptic_disc.num_subdomains)]
initial_basis = elliptic_disc.solution_space.empty()._blocks
for basis in (('1',) if config['initial_basis'] == 0 else ('1', 'x[0]', 'x[1]', 'x[0]*x[1]')):
    vec = make_listvectorarray(wrapper[example.project(basis)])
    vec = BlockVectorArray([elliptic_disc.localize_vector(vec, ss)
                            for ss in np.arange(elliptic_disc.num_subdomains)])
    initial_basis, _ = gram_schmidt_block_basis_extension(initial_basis, vec, product=local_products)
f_h = elliptic_disc.rhs.as_vector(1.)
f_h = elliptic_disc.l2_product.apply_inverse(f_h)
initial_basis, _ = gram_schmidt_block_basis_extension(initial_basis, f_h, product=local_products)
detailed_data['initial_basis'] = initial_basis

greedy_data = reduce_pod_greedy(config, detailed_data, training_samples)

add_values(greedy_max_err_mus=greedy_data['max_err_mus'],
           greedy_max_errs=greedy_data['max_errs'],
           greedy_extensions=greedy_data['extensions'],
           greedy_basis_sizes=[len(local_RB) for local_RB in greedy_data['basis']],
           training_samples=training_samples)
print(' ')

logger.info('creating reference discretization ...')
config_ref = dict(config)
config_ref['dune_num_elements'] = '[200 40]'
config_ref['dune_num_subdomains'] = '[1 1]'
config_ref['nt'] = 2*config['nt']
reference_data = prepare(config_ref)
parabolic_disc_ref, prolongator, bochner_norms_ref = (
        reference_data['parabolic_disc'],
        reference_data['prolongator'],
        reference_data['bochner_norms'])
logger.info('  discretization has {} DoFs'.format(parabolic_disc_ref.solution_space.dim))
print(' ')

logger.info('error analysis ...')
rds, rcs = greedy_data['rds'], greedy_data['rcs']
test_samples = list(parabolic_disc.parameter_space.sample_randomly(config['num_test_samples']))
add_values(test_samples=test_samples,
           num_DoFs=[rd.solution_space.dim for rd in rds])
for mu in test_samples:
    parabolic_disc_ref.solve(mu)
reference_error_computer = DetailedAgainstReference(parabolic_disc_ref,
                                                    prolongator,
                                                    elliptic_disc,
                                                    partial(bochner_norms_ref['elliptic_penalty'], mu=detailed_data['mu_bar']))
max_ests = []
max_errs = []
for rd, rc in zip(rds, rcs):
    ests = []
    errs = []
    for mu in test_samples:
        U_red = rd.solve(mu)
        ests.append(rd.estimate(U_red, mu=mu))
        U_h = rc.reconstruct(U_red)
        errs.append(reference_error_computer.estimate(U_h, mu, parabolic_disc))
    max_ests.append(np.max(ests))
    max_errs.append(np.max(errs))

add_values(test_max_errs=max_errs,
           test_max_ests=max_ests,
           test_num_DoFs=[rd.solution_space.dim for rd in rds])

