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
    #module.add_container('std::vector< int >', 'int', 'vector')
    module.add_container('std::vector< double >', 'double', 'list')
    #module.add_function('id', pybindgen.retval('std::string'), [])
    module.add_function('writeDescriptionFile', None, [])
    #module.add_function('run', pybindgen.retval('int'), [])
    # DuneVector
    DuneVector = module.add_class('DuneVector')
    DuneVector.add_constructor([])
    DuneVector.add_constructor([pybindgen.param('int', 'ss')])
    DuneVector.add_constructor([pybindgen.param('DuneVector *', 'other', transfer_ownership=False)])
    #ColMajorDenseMatrix.add_binary_numeric_operator('+')
    #ColMajorDenseMatrix.add_method('data', pybindgen.retval('double *', caller_owns_return=True), [], is_const=True)
    #ColMajorDenseMatrix.add_method('shape', pybindgen.retval('std::vector< int >'), [], is_const=True)
    DuneVector.add_method('len', pybindgen.retval('int'), [], is_const=True)
    DuneVector.add_method('dot', pybindgen.retval('double'), [pybindgen.param('DuneVector *', 'other', transfer_ownership=False)], is_const=True)
    DuneVector.add_method('scale', None, [pybindgen.param('double', 'scale')], is_const=True)
    DuneVector.add_method('add', pybindgen.retval('DuneVector *', caller_owns_return=True), [pybindgen.param('DuneVector *', 'other', transfer_ownership=False)], is_const=True)
    # DuneOperator
    DuneOperator = module.add_class('DuneOperator')
    DuneOperator.add_method('apply',
                        pybindgen.retval('DuneVector *', caller_owns_return=True),
                        [pybindgen.param('DuneVector *', 'vector', transfer_ownership=False)],
                        is_const=True)
    DuneOperator.add_method('apply2',
                        pybindgen.retval('double'),
                        [pybindgen.param('DuneVector *', 'vectorOne', transfer_ownership=False),
                         pybindgen.param('DuneVector *', 'vectorTwo', transfer_ownership=False)],
                        is_const=True)
    module.add_container('std::vector< DuneOperator * >', pybindgen.retval('DuneOperator *', caller_owns_return=True), 'list')
    # LinearEllipticExampleCG
    LinearEllipticExampleCG = module.add_class('LinearEllipticExampleCG')
    #LinearEllipticExampleCG.add_constructor([pybindgen.param('std::vector< std::string >', 'arg')])
    LinearEllipticExampleCG.add_constructor([])
    LinearEllipticExampleCG.add_method('paramSize', pybindgen.retval('int'), [], is_const=True)
    LinearEllipticExampleCG.add_method('solve',
                                       pybindgen.retval('DuneVector *', caller_owns_return=True),
                                       [pybindgen.param('std::vector< double >', 'mu')],
                                       is_const=True)
    LinearEllipticExampleCG.add_method('operators',
                                       pybindgen.retval('std::vector< DuneOperator * >'),
                                       [],
                                       is_const=True)

    with open(generator_filename, 'wb') as output:
        module.generate(FileCodeSink(output))

    ## check if generate() was successfull
    #print os.stat(generator_fn).st_size

if __name__ == '__main__':
    generate_my_module(sys.argv[1], sys.argv[2], sys.argv[3:])
