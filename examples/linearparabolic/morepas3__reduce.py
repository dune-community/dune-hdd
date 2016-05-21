#!/usr/bin/env python2
#
# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import division, print_function

from scipy.sparse import coo_matrix, bmat

from pymor.algorithms.basisextension import pod_basis_extension
from pymor.algorithms.greedy import greedy
from pymor.operators.constructions import LincombOperator
from pymor.operators.numpy import NumpyMatrixOperator
from pymor.playground.algorithms.blockbasisextension import pod_block_basis_extension
from pymor.playground.operators.block import BlockOperator
from pymor.playground.reductors import GenericBlockRBReconstructor
from pymor.reductors.basic import reduce_generic_rb
from pymor.vectorarrays.block import BlockVectorArray
from pymor.vectorarrays.numpy import NumpyVectorArray


class Reconstructor(object):

    def __init__(self, disc, RB):
        self._disc = disc
        self._RB = RB

    def reconstruct(self, U):
        if U.dim == 0:
            return self._disc.globalize_vectors(self._disc.solution_space.zeros(len(U)))
        else:
            assert U.dim == np.sum(len(RB) for RB in self._RB)
            local_lens = [len(RB) for RB in self._RB]
            local_starts = np.insert(np.cumsum(local_lens), 0, 0)[:-1]
            localized_U = [NumpyVectorArray(U._array[:, local_starts[ii]:(local_starts[ii] + local_lens[ii])])
                           for ii in np.arange(len(self._RB))]
            U = GenericBlockRBReconstructor(self._RB).reconstruct(BlockVectorArray(localized_U))
            return self._disc.globalize_vectors(U)

    def restricted_to_subbasis(self, dim):
        if not isinstance(dim, tuple):
            dim = len(self.RB)*[dim]
        assert all([dd <= len(rb) for dd, rb in izip(dim, self.RB)])
        return Reconstructor(self._disc, [rb.copy(ind=range(dd)) for rb, dd in izip(self.RB, dim)])


def reductor(discretization, RB, vector_product=None, disable_caching=True, extends=None):
    if RB is None:
        RB = [elliptic_LRBMS_disc.local_operator(ss).source.empty() for ss in np.arange(elliptic_LRBMS_disc.num_subdomains)]
    rd, rc, reduction_data = reduce_generic_rb(elliptic_LRBMS_disc, RB, vector_product, disable_caching, extends)

    def unblock_op(op):
        assert op._blocks[0][0] is not None
        if isinstance(op._blocks[0][0], LincombOperator):
            coefficients = op._blocks[0][0].coefficients
            operators = [None for kk in np.arange(len(op._blocks[0][0].operators))]
            for kk in np.arange(len(op._blocks[0][0].operators)):
                ops = [[op._blocks[ii][jj].operators[kk]
                        if op._blocks[ii][jj] is not None else None
                        for jj in np.arange(op.num_source_blocks)]
                       for ii in np.arange(op.num_range_blocks)]
                operators[kk] = unblock_op(BlockOperator(ops))
            return LincombOperator(operators=operators, coefficients=coefficients)
        else:
            assert all(all([isinstance(block, NumpyMatrixOperator) if block is not None else True
                           for block in row])
                       for row in op._blocks)
            if op.source.dim == 0 and op.range.dim == 0:
                return NumpyMatrixOperator(np.zeros((0, 0)))
            elif op.source.dim == 1:
                mat = np.concatenate([op._blocks[ii][0]._matrix
                                      for ii in np.arange(op.num_range_blocks)],
                                     axis=1)
            elif op.range.dim == 1:
                mat = np.concatenate([op._blocks[0][jj]._matrix
                                      for jj in np.arange(op.num_source_blocks)],
                                     axis=1)
            else:
                mat = bmat([[coo_matrix(op._blocks[ii][jj]._matrix)
                             if op._blocks[ii][jj] is not None else coo_matrix((op._range_dims[ii], op._source_dims[jj]))
                             for jj in np.arange(op.num_source_blocks)]
                            for ii in np.arange(op.num_range_blocks)]).toarray()
            return NumpyMatrixOperator(mat)

    class DetailedEstimator(object):

        def estimate(self, U, mu, discretization):
            U_rec = Reconstructor(elliptic_LRBMS_disc, RB).reconstruct(U)
            U_prol = prolong(reference_example, elliptic_LRBMS_disc._impl, U_rec)
            return reference_norm(parabolic_reference_disc.solve(mu) - U_prol)

    reduced_op = unblock_op(rd.operator)
    reduced_rhs = unblock_op(rd.rhs)
    return (InstationaryDiscretization(T=config['end_time'],
                                       initial_data=reduced_op.source.zeros(1),
                                       operator=reduced_op,
                                       rhs=unblock_op(rd.rhs),
                                       mass=unblock_op(rd.products['l2']),
                                       time_stepper=ImplicitEulerTimeStepper(config['nt']),
                                       products={kk: unblock_op(rd.products[kk]) for kk in rd.products.keys()},
                                       operators={kk: unblock_op(rd.operators[kk])
                                                  for kk in rd.operators.keys() if kk != 'operator'},
                                       functionals={kk: unblock_op(rd.functionals[kk])
                                                    for kk in rd.functionals.keys() if kk != 'rhs'},
                                       vector_operators={kk: unblock_op(rd.vector_operators[kk])
                                                         for kk in rd.vector_operators.keys()},
                                       parameter_space=rd.parameter_space,
                                       estimator=DetailedEstimator(),
                                       cache_region='disk',
                                       name='reduced discretization'),
            Reconstructor(elliptic_LRBMS_disc, RB),
            reduction_data)


def extension(basis, U):
    if not isinstance(U, BlockVectorArray):
        U = BlockVectorArray([elliptic_LRBMS_disc.localize_vector(U, ii)
                              for ii in np.arange(elliptic_LRBMS_disc.num_subdomains)])
    return pod_block_basis_extension(basis,
                                     U,
                                     count=1,
                                     product=[elliptic_LRBMS_disc.local_product(ss, config['extension_product'])
                                              for ss in np.arange(elliptic_LRBMS_disc.num_subdomains)])


def reduce():
    greedy_data = greedy(parabolic_disc,
                         reductor,
                         training_samples,
                         use_estimator=True,
                         error_norm=None,
                         extension_algorithm=extension,
                         max_extensions=config['max_rb_size'],
                         target_error=config['target_error'])
    rd, rc = greedy_data['reduced_discretization'], greedy_data['reconstructor']
    RB = greedy_data['basis']
