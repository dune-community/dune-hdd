#! /usr/bin/env python

import sys, os
import pybindgen
from pybindgen import FileCodeSink
#from pybindgen.typehandlers import base as typehandlers
#from pybindgen import ReturnValue, Parameter, Module, Function, FileCodeSink
#from pybindgen import CppMethod, CppConstructor, CppClass, Enum

def generate_my_module(inputdir, outputdir, includedirs):
    includedirs = includedirs[0].split(';')
    generator_filename = os.path.join(outputdir, 'dunelinearellipticcg2dsgrid.cc')
    module = pybindgen.Module('dunelinearellipticcg2dsgrid')
    module.add_include('"detailed_discretizations.hh"')
    #module.add_function('run', None, [])
    #namespace = module.add_cpp_namespace('Dune').add_cpp_namespace('Stuff').add_cpp_namespace('LA').add_cpp_namespace('Container')
    module.add_container('std::vector< std::string >', 'std::string', 'list')
    module.add_container('std::vector< int >', 'int', 'vector')
    module.add_container('std::vector< double >', 'double', 'vector')
    module.add_function('writeDescriptionFile', None, [pybindgen.param('std::string', 'filename'),
                                                       pybindgen.param('std::string', '_id')])
    # ColMajorDenseMatrix
    ColMajorDenseMatrix = module.add_class('ColMajorDenseMatrix')
    ColMajorDenseMatrix.add_constructor([])
    ColMajorDenseMatrix.add_constructor([pybindgen.param('int', '_rows'), pybindgen.param('int', '_cols')])
    ColMajorDenseMatrix.add_constructor([pybindgen.param('ColMajorDenseMatrix', 'other')])
    ColMajorDenseMatrix.add_binary_numeric_operator('+')
    ColMajorDenseMatrix.add_method('data', pybindgen.retval('double *', caller_owns_return=True), [], is_const=True)
    ColMajorDenseMatrix.add_method('shape', pybindgen.retval('std::vector< int >'), [], is_const=True)
    ColMajorDenseMatrix.add_method('len', pybindgen.retval('int'), [], is_const=True)
    # Operator
    Operator = module.add_class('Operator')
    Operator.add_method('apply',
                        pybindgen.retval('ColMajorDenseMatrix'),
                        [pybindgen.param('ColMajorDenseMatrix', 'vectors')],
                        is_const=True)
    Operator.add_method('apply2',
                        pybindgen.retval('ColMajorDenseMatrix'),
                        [pybindgen.param('ColMajorDenseMatrix', 'vectors'),
                         pybindgen.param('ColMajorDenseMatrix', 'vectors')],
                        is_const=True)
    # LinearEllipticExampleCG
    LinearEllipticExampleCG = module.add_class('LinearEllipticExampleCG')
    LinearEllipticExampleCG.add_constructor([pybindgen.param('std::vector< std::string >', 'arg')])
    LinearEllipticExampleCG.add_method('parametric', pybindgen.retval('bool'), [], is_const=True)
    LinearEllipticExampleCG.add_method('solve',
                                       pybindgen.retval('ColMajorDenseMatrix *', caller_owns_return=True),
                                       [])
    LinearEllipticExampleCG.add_method('getOperator',
                                       pybindgen.retval('Operator *', caller_owns_return=True),
                                       [])

    with open(generator_filename, 'wb') as output:
        module.generate(FileCodeSink(output))

    ## check if generate() was successfull
    #print os.stat(generator_fn).st_size

if __name__ == '__main__':
    generate_my_module(sys.argv[1], sys.argv[2], sys.argv[3:])
