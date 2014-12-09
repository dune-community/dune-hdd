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
    ssize_t = CONFIG_H['DUNE_STUFF_SSIZE_T']
    '''injects the user code into the module'''
    # first the discretization
    GridType = 'Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming >'
    #GridType = 'Dune::SGrid< 2, 2 >'
    RangeFieldType = 'double'
    dimRange = '1'
    polOrder = '1'
    la_backend = 'Dune::Stuff::LA::ChooseBackend::eigen_sparse'
    if la_backend == 'Dune::Stuff::LA::ChooseBackend::eigen_sparse':
        MatrixType = 'Dune::Stuff::LA::EigenRowMajorSparseMatrix< ' + RangeFieldType + ' >'
        VectorType = 'Dune::Stuff::LA::EigenDenseVector< ' + RangeFieldType + ' >'
    elif la_backend == 'Dune::Stuff::LA::ChooseBackend::istl_sparse':
        MatrixType = 'Dune::Stuff::LA::IstlRowMajorSparseMatrix< ' + RangeFieldType + ' >'
        VectorType = 'Dune::Stuff::LA::IstlDenseVector< ' + RangeFieldType + ' >'
    space_backend = 'Dune::GDT::ChooseSpaceBackend::fem'
    OperatorType = 'Dune::Pymor::Operators::LinearAffinelyDecomposedContainerBased< ' + MatrixType + ', ' + VectorType + ' >'
    ProductType = OperatorType
    FunctionalType = 'Dune::Pymor::Functionals::LinearAffinelyDecomposedVectorBased< ' + VectorType + ' >'
    DiscretizationName = 'Dune::HDD::LinearElliptic::Discretizations::BlockSWIPDG'
    DiscretizationFullName = (DiscretizationName + '< '
                              + GridType + ', '
                              + RangeFieldType + ', '
                              + dimRange + ', ' + polOrder + ', '
                              + la_backend + '>')
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
    # discretization = inject_StationaryDiscretizationImplementation(
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
    # then add the example
    def add_example(name):
        Example = module.add_class(name, template_parameters=[GridType], custom_name=name)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors'),
                                 param('const std::string', 'info_color'),
                                 param('const std::string', 'debug_color'),
                                 param('const std::string', 'warn_color')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors'),
                                 param('const std::string', 'info_color'),
                                 param('const std::string', 'debug_color')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors'),
                                 param('const std::string', 'info_color')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings'),
                                 param('const bool', 'enable_colors')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels'),
                                 param('const bool', 'enable_warnings')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products'),
                                 param('const ' + ssize_t, 'info_log_levels'),
                                 param('const ' + ssize_t, 'debug_log_levels')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products'),
                                 param('const ' + ssize_t, 'info_log_levels')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers'),
                                 param('const std::vector< std::string >', 'products')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements'),
                                 param('const ' + ssize_t, 'oversampling_layers')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning'),
                                 param('const ' + ssize_t, 'num_refinements')],
                                throw=exceptions)
        Example.add_constructor([param('const std::string', 'partitioning')],
                                throw=exceptions)
        Example.add_constructor([], throw=exceptions)
        Example.add_method('discretization_and_return_ptr',
                           retval(DiscretizationFullName + ' *', caller_owns_return=True),
                           [], is_const=True, throw=exceptions,
                           custom_name='discretization')
        Example.add_method('visualize',
                           None,
                           [param('const std::string&', 'filename_prefix')],
                           is_const=True, throw=exceptions)
        Example.add_method('project',
                           retval(VectorType),
                           [param('const std::string', 'expression')], is_const=True, throw=exceptions)
        Example.add_method('compute_error',
                           retval(RangeFieldType),
                           [param('const ' + VectorType + '&', 'solution'),
                            param('const std::string', 'product_type'),
                            param('const Dune::Pymor::Parameter', 'mu'),
                            param('const Dune::Pymor::Parameter', 'mu_product')],
                           is_const=True, throw=exceptions)
        Example.add_method('compute_error',
                           retval(RangeFieldType),
                           [param('const ' + VectorType + '&', 'solution'),
                            param('const std::string', 'product_type'),
                            param('const Dune::Pymor::Parameter', 'mu')],
                           is_const=True, throw=exceptions)
        Example.add_method('compute_error',
                           retval(RangeFieldType),
                           [param('const ' + VectorType + '&', 'solution'),
                            param('const std::string', 'product_type')],
                           is_const=True, throw=exceptions)
        Example.add_method('compute_error',
                           retval(RangeFieldType),
                           [param('const ' + VectorType + '&', 'solution')],
                           is_const=True, throw=exceptions)
        Example.add_method('available_estimators',
                           retval('std::vector< std::string >'),
                           [],
                           is_const=True, throw=exceptions)
        Example.add_method('estimate',
                           retval(RangeFieldType),
                           [param('const ' + VectorType + '&', 'vector'),
                            param('const std::string', 'type'),
                            param('const Dune::Pymor::Parameter', 'mu_hat'),
                            param('const Dune::Pymor::Parameter', 'mu_bar'),
                            param('const Dune::Pymor::Parameter', 'mu')],
                           is_const=True, throw=exceptions)
        Example.add_method('estimate',
                           retval(RangeFieldType),
                           [param('const ' + VectorType + '&', 'vector'),
                            param('const std::string', 'type')],
                           is_const=True, throw=exceptions)
        Example.add_method('available_local_estimators',
                           retval('std::vector< std::string >'),
                           [],
                           is_const=True, throw=exceptions)
        Example.add_method('estimate_local',
                           retval('std::vector< ' + RangeFieldType + ' >'),
                           [param('const ' + VectorType + '&', 'vector'),
                            param('const std::string', 'type'),
                            param('const Dune::Pymor::Parameter', 'mu_hat'),
                            param('const Dune::Pymor::Parameter', 'mu_bar'),
                            param('const Dune::Pymor::Parameter', 'mu')],
                           is_const=True, throw=exceptions)
        Example.add_method('estimate_local',
                           retval('std::vector< ' + RangeFieldType + ' >'),
                           [param('const ' + VectorType + '&', 'vector'),
                            param('const std::string', 'type')],
                           is_const=True, throw=exceptions)
        Example.add_method('pb_project_global_to_oversampled',
                           retval(VectorType + '*', caller_owns_return=True),
                           [param('const ' + VectorType + '&', 'global_vector'),
                            param('const ' + ssize_t, 'subdomain')],
                           is_const=True, throw=exceptions,
                           custom_name='project_global_to_oversampled')
        Example.add_method('pb_project_global_to_local',
                           retval(VectorType + '*', caller_owns_return=True),
                           [param('const ' + VectorType + '&', 'global_vector'),
                            param('const ' + ssize_t, 'subdomain')],
                           is_const=True, throw=exceptions,
                           custom_name='project_global_to_local')
        Example.add_method('pb_project_oversampled_to_local',
                           retval(VectorType + '*', caller_owns_return=True),
                           [param('const ' + VectorType + '&', 'oversampled_vector'),
                            param('const ' + ssize_t, 'subdomain')],
                           is_const=True, throw=exceptions,
                           custom_name='project_oversampled_to_local')
        Example.add_method('solve_for_local_correction',
                           retval(VectorType),
                           [param('std::vector< ' + VectorType + ' >&', 'local_vectors'),
                            param('const ' + ssize_t, 'subdomain'),
                            param('const Dune::Pymor::Parameter', 'mu')],
                           is_const=True, throw=exceptions)
        Example.add_method('solve_for_local_correction',
                           retval(VectorType),
                           [param('std::vector< ' + VectorType + ' >&', 'local_vectors'),
                            param('const ' + CONFIG_H['DUNE_STUFF_SSIZE_T'], 'subdomain')],
                           is_const=True, throw=exceptions)
        Example.add_method('solve_oversampled',
                           retval(VectorType),
                           [param('const ' + ssize_t, 'subdomain'),
                            param('const std::string&', 'boundary_value_type'),
                            param('const ' + VectorType + '&', 'boundary_value_type'),
                            param('const Dune::Pymor::Parameter', 'mu')],
                           is_const=True, throw=exceptions)
        Example.add_method('solve_oversampled',
                           retval(VectorType),
                           [param('const ' + ssize_t, 'subdomain'),
                            param('const std::string&', 'boundary_value_type'),
                            param('const ' + VectorType + '&', 'boundary_value_type')],
                           is_const=True, throw=exceptions)
        Example.add_method('alpha',
                           retval(RangeFieldType),
                           [param('const Dune::Pymor::Parameter&', 'mu_1'),
                            param('const Dune::Pymor::Parameter&', 'mu_2')],
                           is_const=True, throw=exceptions)
        Example.add_method('gamma',
                           retval(RangeFieldType),
                           [param('const Dune::Pymor::Parameter&', 'mu_1'),
                            param('const Dune::Pymor::Parameter&', 'mu_2')],
                           is_const=True, throw=exceptions)
    add_example('Spe10Model1Example')
    add_example('OS2014Example')
    add_example('FiveSpotExample')
    add_example('LocalThermalblockExample')


if __name__ == '__main__':
    # prepare the module
    module, pybindgen_filename, config_h_filename = prepare_python_bindings(sys.argv[1:])
    # add all of libdunepymor
    module, exceptions, interfaces, CONFIG_H = inject_lib_dune_pymor(module, config_h_filename)
    # add example user code (see above)
    inject_Example(module, exceptions, interfaces, CONFIG_H)
    # and finally write the pybindgen .cc file
    finalize_python_bindings(module, pybindgen_filename)
