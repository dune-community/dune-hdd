#! /usr/bin/env python
# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import sys
from pybindgen import param, retval

from dune.pymor.core import prepare_python_bindings, inject_lib_dune_pymor, finalize_python_bindings
from dune.pymor.discretizations import inject_StationaryDiscretizationImplementation
from dune.pymor.discretizations import inject_StationaryMultiscaleDiscretizationImplementation


def inject_Example(module, exceptions, interfaces, CONFIG_H):
    '''injects the user code into the module'''
    ssize_t = CONFIG_H['DUNE_STUFF_SSIZE_T']
    HAVE_DUNE_PDELAB = CONFIG_H['HAVE_DUNE_PDELAB']
    HAVE_DUNE_FEM    = CONFIG_H['HAVE_DUNE_FEM']
    HAVE_DUNE_ISTL   = CONFIG_H['HAVE_DUNE_ISTL']
    HAVE_EIGEN       = CONFIG_H['HAVE_EIGEN']
    HAVE_ALUGRID     = CONFIG_H['HAVE_ALUGRID']
    HAVE_DUNE_SPGRID = CONFIG_H['HAVE_DUNE_SPGRID']
    RangeFieldType = 'double'
    def add_example(GridType, space_backend, la_backend, name):
        # build all types needed for the discretization
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
        LocalDiscretizationName = 'Dune::HDD::LinearElliptic::Discretizations::SWIPDG'
        local_layer = 'Dune::Stuff::Grid::ChooseLayer::local'
        LocalDiscretizationFullName = (LocalDiscretizationName
                                       + '< ' + GridType + ', '
                                       + local_layer + ', '
                                       + RangeFieldType + ', '
                                       + dimRange + ', ' + polOrder + ', '
                                       + space_backend + ', '
                                       + la_backend + ' >')
        inject_StationaryDiscretizationImplementation(module, exceptions, interfaces, CONFIG_H,
                                                      LocalDiscretizationName,
                                                      Traits={'VectorType': VectorType,
                                                              'OperatorType': OperatorType,
                                                              'FunctionalType': FunctionalType,
                                                              'ProductType': ProductType},
                                                      template_parameters=[GridType, local_layer, RangeFieldType,
                                                                           dimRange, polOrder, space_backend, la_backend])
        oversampled_layer = 'Dune::Stuff::Grid::ChooseLayer::local_oversampled'
        OversampledDiscretizationFullName = (LocalDiscretizationName
                                             + '< ' + GridType + ', '
                                             + oversampled_layer + ', '
                                             + RangeFieldType + ', '
                                             + dimRange + ', ' + polOrder + ', '
                                             + space_backend + ', '
                                             + la_backend + ' >')
        inject_StationaryDiscretizationImplementation(module, exceptions, interfaces, CONFIG_H,
                                                      LocalDiscretizationName,
                                                      Traits={'VectorType': VectorType,
                                                              'OperatorType': OperatorType,
                                                              'FunctionalType': FunctionalType,
                                                              'ProductType': ProductType},
                                                      template_parameters=[GridType, oversampled_layer, RangeFieldType,
                                                                           dimRange, polOrder, space_backend, la_backend])
        DiscretizationName = 'Dune::HDD::LinearElliptic::Discretizations::BlockSWIPDG'
        DiscretizationType = (DiscretizationName + '< '
                              + GridType + ', '
                              # + GridType + ', ' + grid_layer + ', '
                              + RangeFieldType + ', '
                              + dimRange + ', ' + polOrder + ', '
                              + la_backend + '>')
        inject_StationaryMultiscaleDiscretizationImplementation(
                module, exceptions, interfaces, CONFIG_H,
                DiscretizationName,
                Traits={'VectorType': VectorType,
                        'OperatorType': OperatorType,
                        'FunctionalType': FunctionalType,
                        'ProductType': ProductType,
                        'LocalDiscretizationType': LocalDiscretizationFullName,
                        'OversampledDiscretizationType': OversampledDiscretizationFullName},
                template_parameters=[GridType, RangeFieldType,
                                     dimRange, polOrder, la_backend])
        # then create the example
        Example = module.add_class('PbGenericLinearellipticMultiscaleExample', template_parameters=[GridType, space_backend, la_backend], custom_name=name)
        Example.add_method('logger_options',
                           retval('Dune::Stuff::Common::Configuration'),
                           [], is_static=True, throw=exceptions)
        Example.add_method('grid_options',
                           retval('std::vector< std::string >'),
                           [], is_static=True, throw=exceptions)
        Example.add_method('grid_options',
                           retval('Dune::Stuff::Common::Configuration'),
                           [param('const std::string&', 'type')], is_static=True, throw=exceptions)
        Example.add_method('boundary_options',
                           retval('std::vector< std::string >'),
                           [], is_static=True, throw=exceptions)
        Example.add_method('boundary_options',
                           retval('Dune::Stuff::Common::Configuration'),
                           [param('const std::string&', 'type')], is_static=True, throw=exceptions)
        Example.add_method('problem_options',
                           retval('std::vector< std::string >'),
                           [], is_static=True, throw=exceptions)
        Example.add_method('problem_options',
                           retval('Dune::Stuff::Common::Configuration'),
                           [param('const std::string&', 'type')], is_static=True, throw=exceptions)
        Example.add_method('solver_options',
                           retval('std::vector< std::string >'),
                           [], is_static=True, throw=exceptions)
        Example.add_method('solver_options',
                           retval('Dune::Stuff::Common::Configuration'),
                           [param('const std::string&', 'type')], is_static=True, throw=exceptions)
        Example.add_constructor([param('const Dune::Stuff::Common::Configuration&', 'logger_cfg'),
                                 param('const Dune::Stuff::Common::Configuration&', 'grid_cfg'),
                                 param('const Dune::Stuff::Common::Configuration&', 'boundary_cfg'),
                                 param('const Dune::Stuff::Common::Configuration&', 'problem_cfg')],
                                throw=exceptions)
        Example.add_method('pb_discretization_and_return_ptr',
                           retval(DiscretizationType + ' *', caller_owns_return=True),
                           [], is_const=True, throw=exceptions,
                           custom_name='discretization')
        Example.add_method('project',
                           retval(VectorType),
                           [param('const std::string', 'expression')], is_const=True, throw=exceptions)
        Example.add_method('prolong',
                           retval(VectorType),
                           [param('const ' + DiscretizationType + '&', 'source_disc'),
                            param('const ' + VectorType + '&', 'source_vec')],
                           is_const=True, throw=exceptions)
        Example.add_method('oswald_interpolate',
                           retval(VectorType),
                           [param('const ' + VectorType + '&', 'vector')],
                           is_const=True, throw=exceptions)
        Example.add_method('alpha',
                           retval('double'),
                           [param('const Dune::Pymor::Parameter&', 'mu_1'),
                            param('const Dune::Pymor::Parameter&', 'mu_2')],
                           is_const=True, throw=exceptions)
        Example.add_method('gamma',
                           retval('double'),
                           [param('const Dune::Pymor::Parameter&', 'mu_1'),
                            param('const Dune::Pymor::Parameter&', 'mu_2')],
                           is_const=True, throw=exceptions)
        Example.add_method('min_diffusion_ev',
                           retval('double'),
                           [param('const Dune::Pymor::Parameter&', 'mu')],
                           is_const=True, throw=exceptions)
        Example.add_method('max_diffusion_ev',
                           retval('double'),
                           [param('const Dune::Pymor::Parameter&', 'mu')],
                           is_const=True, throw=exceptions)
        Example.add_method('visualize',
                           None,
                           [param('const std::string&', 'filename_prefix')],
                           is_const=True, throw=exceptions)
    SGrid1d = 'Dune::SGrid< 1, 1 >'
    SGrid2d = 'Dune::SGrid< 2, 2 >'
    SGrid3d = 'Dune::SGrid< 3, 3 >'
    if HAVE_ALUGRID:
        AluGridConform2d = 'Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming >'
        AluGridConform3d = 'Dune::ALUGrid< 3, 3, Dune::simplex, Dune::conforming >'
    la_backend_eigen = 'Dune::Stuff::LA::ChooseBackend::eigen_sparse'
    la_backend_istl  = 'Dune::Stuff::LA::ChooseBackend::istl_sparse'
    space_backend_pdelab = 'Dune::GDT::ChooseSpaceBackend::pdelab'
    space_backend_fem    = 'Dune::GDT::ChooseSpaceBackend::fem'
    if HAVE_DUNE_FEM and HAVE_DUNE_ISTL:
        add_example(SGrid1d, space_backend_fem, la_backend_istl, 'GenericLinearellipticMultiscaleExample_1dSGrid_fem_istl')
        add_example(SGrid2d, space_backend_fem, la_backend_istl, 'GenericLinearellipticMultiscaleExample_2dSGrid_fem_istl')
        add_example(SGrid3d, space_backend_fem, la_backend_istl, 'GenericLinearellipticMultiscaleExample_3dSGrid_fem_istl')
        if HAVE_ALUGRID:
            add_example(AluGridConform2d, space_backend_fem, la_backend_istl,
                        'GenericLinearellipticMultiscaleExample_2dAluConformGrid_fem_istl')
            add_example(AluGridConform3d, space_backend_fem, la_backend_istl,
                        'GenericLinearellipticMultiscaleExample_3dAluConformGrid_fem_istl')
    if HAVE_DUNE_FEM and HAVE_EIGEN:
        add_example(SGrid1d, space_backend_fem, la_backend_eigen, 'GenericLinearellipticMultiscaleExample_1dSGrid_fem_eigen')
        add_example(SGrid2d, space_backend_fem, la_backend_eigen, 'GenericLinearellipticMultiscaleExample_2dSGrid_fem_eigen')
        add_example(SGrid3d, space_backend_fem, la_backend_eigen, 'GenericLinearellipticMultiscaleExample_3dSGrid_fem_eigen')
        if HAVE_ALUGRID:
            add_example(AluGridConform2d, space_backend_fem, la_backend_eigen,
                        'GenericLinearellipticMultiscaleExample_2dAluConformGrid_fem_eigen')
            add_example(AluGridConform3d, space_backend_fem, la_backend_eigen,
                        'GenericLinearellipticMultiscaleExample_3dAluConformGrid_fem_eigen')


if __name__ == '__main__':
    # prepare the module
    module, pybindgen_filename, config_h_filename = prepare_python_bindings(sys.argv[1:])
    # add all of libdunepymor
    module, exceptions, interfaces, CONFIG_H = inject_lib_dune_pymor(module, config_h_filename)
    # add example user code (see above)
    inject_Example(module, exceptions, interfaces, CONFIG_H)
    # and finally write the pybindgen .cc file
    finalize_python_bindings(module, pybindgen_filename)
