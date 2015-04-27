#! /usr/bin/env python
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import sys
from pybindgen import param, retval

from dune.pymor.core import prepare_python_bindings, inject_lib_dune_pymor, finalize_python_bindings
from dune.pymor.discretizations import inject_StationaryDiscretizationImplementation


def inject_Example(module, exceptions, interfaces, CONFIG_H):
    '''injects the user code into the module'''
    ssize_t = CONFIG_H['DUNE_STUFF_SSIZE_T']
    HAVE_DUNE_PDELAB = CONFIG_H['HAVE_DUNE_PDELAB']
    HAVE_DUNE_FEM    = CONFIG_H['HAVE_DUNE_FEM']
    HAVE_DUNE_ISTL   = CONFIG_H['HAVE_DUNE_ISTL']
    HAVE_EIGEN       = CONFIG_H['HAVE_EIGEN']
    def add_example(GridType, space_backend, la_backend, name):
        # build all types needed for the discretization
        RangeFieldType = 'double'
        dimRange = '1'
        polOrder = '1'
        grid_layer = 'Dune::Stuff::Grid::ChooseLayer::leaf'
        MatrixType = 'Dune::Stuff::LA::'
        VectorType = 'Dune::Stuff::LA::'
        if 'eigen_sparse' in la_backend:
            MatrixType += 'EigenRowMajorSparseMatrix'
            VectorType += 'EigenDenseVector'
        elif 'istl_sparse' in la_backend:
            MatrixType += 'IstlRowMajorSparseMatrix'
            VectorType += 'IstlDenseVector'
        MatrixType += '< ' + RangeFieldType + ' >'
        VectorType += '< ' + RangeFieldType + ' >'
        OperatorType = 'Dune::Pymor::Operators::LinearAffinelyDecomposedContainerBased< ' + MatrixType + ', ' + VectorType + ' >'
        ProductType = OperatorType
        FunctionalType = 'Dune::Pymor::Functionals::LinearAffinelyDecomposedVectorBased< ' + VectorType + ' >'
        DiscretizationName = 'Dune::HDD::LinearElliptic::Discretizations::CG'
        DiscretizationType = (DiscretizationName + '< '
                              + GridType + ', ' + grid_layer + ', '
                              + RangeFieldType + ', '
                              + dimRange + ', ' + polOrder + ', '
                              + space_backend + ', ' + la_backend + '>')
        inject_StationaryDiscretizationImplementation(module, exceptions, interfaces, CONFIG_H,
                                                      DiscretizationName,
                                                      Traits={'VectorType': VectorType,
                                                              'OperatorType': OperatorType,
                                                              'FunctionalType': FunctionalType,
                                                              'ProductType': ProductType},
                                                      template_parameters=[GridType, grid_layer, RangeFieldType,
                                                                           dimRange, polOrder, space_backend, la_backend])
        # then create the example
        Example = module.add_class('PbCgExample', template_parameters=[GridType, space_backend, la_backend], custom_name=name)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors'),
                                 param('const std::string', 'info_color'),
                                 param('const std::string', 'debug_color'),
                                 param('const std::string', 'warn_color')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors'),
                                 param('const std::string', 'info_color'),
                                 param('const std::string', 'debug_color')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors'),
                                 param('const std::string', 'info_color')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements'),
                                 param('const ' + ssize_t, 'info_log_levels')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks'),
                                 param('const std::string&', 'num_grid_elements')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string&', 'num_blocks')],
                                throw=exceptions)
        Example.add_constructor([],
                                throw=exceptions)
        Example.add_method('pb_discretization_and_return_ptr',
                           retval(DiscretizationType + ' *', caller_owns_return=True),
                           [], is_const=True, throw=exceptions,
                           custom_name='discretization')
    YaspGrid2d = 'Dune::YaspGrid< 2 >'
    YaspGrid3d = 'Dune::YaspGrid< 3 >'
    la_backend_eigen = 'Dune::Stuff::LA::ChooseBackend::eigen_sparse'
    la_backend_istl  = 'Dune::Stuff::LA::ChooseBackend::istl_sparse'
    space_backend_pdelab = 'Dune::GDT::ChooseSpaceBackend::pdelab'
    space_backend_fem    = 'Dune::GDT::ChooseSpaceBackend::fem'
    if HAVE_DUNE_PDELAB and HAVE_DUNE_ISTL:
        add_example(YaspGrid2d, space_backend_pdelab, la_backend_istl, 'CG_Thermalblock_2dYaspGrid_pdelab_istl')
        add_example(YaspGrid3d, space_backend_pdelab, la_backend_istl, 'CG_Thermalblock_3dYaspGrid_pdelab_istl')
    if HAVE_DUNE_PDELAB and HAVE_EIGEN:
        add_example(YaspGrid2d, space_backend_pdelab, la_backend_eigen, 'CG_Thermalblock_2dYaspGrid_pdelab_eigen')
        add_example(YaspGrid3d, space_backend_pdelab, la_backend_eigen, 'CG_Thermalblock_3dYaspGrid_pdelab_eigen')
    if HAVE_DUNE_FEM and HAVE_DUNE_ISTL:
        add_example(YaspGrid2d, space_backend_fem, la_backend_istl, 'CG_Thermalblock_2dYaspGrid_fem_istl')
        add_example(YaspGrid3d, space_backend_fem, la_backend_istl, 'CG_Thermalblock_3dYaspGrid_fem_istl')
    if HAVE_DUNE_FEM and HAVE_EIGEN:
        add_example(YaspGrid2d, space_backend_fem, la_backend_eigen, 'CG_Thermalblock_2dYaspGrid_fem_eigen')
        add_example(YaspGrid3d, space_backend_fem, la_backend_eigen, 'CG_Thermalblock_3dYaspGrid_fem_eigen')


if __name__ == '__main__':
    # prepare the module
    module, pybindgen_filename, config_h_filename = prepare_python_bindings(sys.argv[1:])
    # add all of libdunepymor
    module, exceptions, interfaces, CONFIG_H = inject_lib_dune_pymor(module, config_h_filename)
    # add example user code (see above)
    inject_Example(module, exceptions, interfaces, CONFIG_H)
    # and finally write the pybindgen .cc file
    finalize_python_bindings(module, pybindgen_filename)
