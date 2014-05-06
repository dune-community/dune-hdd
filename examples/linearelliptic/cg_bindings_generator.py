#! /usr/bin/env python
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import sys
from pybindgen import param, retval

from dune.pymor.core import prepare_python_bindings, inject_lib_dune_pymor, finalize_python_bindings
from dune.pymor.discretizations import inject_StationaryDiscretizationImplementation


def inject_Example(module, exceptions, interfaces, CONFIG_H):
    '''injects the user code into the module'''
    # first the discretization
    GridType = 'Dune::SGrid< 2, 2 >'
    RangeFieldType = 'double'
    dimRange = '1'
    polOrder = '1'
    MatrixType = 'Dune::Stuff::LA::EigenRowMajorSparseMatrix< ' + RangeFieldType + ' >'
    VectorType = 'Dune::Stuff::LA::EigenDenseVector< ' + RangeFieldType + ' >'
    DiscretizationName = 'Dune::HDD::LinearElliptic::Discretizations::CG'
    DiscretizationFullName = (DiscretizationName + '< '
                              + GridType + ', '
                              + RangeFieldType + ', '
                              + dimRange + ', ' + polOrder + ' >')
    discretization = inject_StationaryDiscretizationImplementation(
        module, exceptions, interfaces, CONFIG_H,
        DiscretizationName,
        Traits={'VectorType': VectorType,
                'OperatorType': 'Dune::Pymor::Operators::LinearAffinelyDecomposedContainerBased< ' + MatrixType + ', ' + VectorType + ' >',
                'FunctionalType': 'Dune::Pymor::Functionals::LinearAffinelyDecomposedVectorBased< ' + VectorType + ' >',
                'ProductType': 'Dune::Pymor::Operators::LinearAffinelyDecomposedContainerBased< ' + MatrixType + ', ' + VectorType + ' > '},
        template_parameters=[GridType, RangeFieldType, dimRange, polOrder])
    # then add the example
    LinearellipticExampleCG = module.add_class('LinearellipticExampleCG',
                                               template_parameters=['Dune::SGrid< 2, 2 >'],
                                               custom_name='LinearellipticExampleCG')
    LinearellipticExampleCG.add_method('static_id',
                                       retval('std::string'),
                                       [], is_const=True, throw=[exceptions['Exception']])
    LinearellipticExampleCG.add_method('write_config_file',
                                       None, [], is_const=True, throw=[exceptions['Exception']])
    LinearellipticExampleCG.add_constructor([], throw=[exceptions['Exception']])
    LinearellipticExampleCG.add_method('initialize', None,
                                       [param('const std::vector< std::string >', 'arguments')],
                                       is_const=True, throw=[exceptions['Exception']])
    LinearellipticExampleCG.add_method('discretization_and_return_ptr',
                                       retval(DiscretizationFullName + ' *', caller_owns_return=True),
                                       [], is_const=True, throw=[exceptions['Exception']],
                                       custom_name='discretization')

if __name__ == '__main__':
    # prepare the module
    module, pybindgen_filename = prepare_python_bindings(sys.argv[1:])
    # add all of libdunepymor
    module, exceptions, interfaces, CONFIG_H = inject_lib_dune_pymor(module)
    # add example user code (see above)
    inject_Example(module, exceptions, interfaces, CONFIG_H)
    # and finally write the pybindgen .cc file
    finalize_python_bindings(module, pybindgen_filename)
