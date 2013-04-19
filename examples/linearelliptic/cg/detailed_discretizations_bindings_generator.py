#! /usr/bin/env python

import sys, os
#, shutil
import pybindgen
#from pybindgen.gccxmlparser import ModuleParser
from pybindgen import FileCodeSink
#from pygccxml.declarations.class_declaration import class_t
#from pybindgen.typehandlers import base as typehandlers
#from pybindgen import ReturnValue, Parameter, Module, Function, FileCodeSink
#from pybindgen import CppMethod, CppConstructor, CppClass, Enum
#from pygccxml.declarations.calldef import free_function_t, member_function_t, constructor_t, calldef_t

def generate_my_module(inputdir, outputdir, includedirs):
    includedirs = includedirs[0].split(';')
    generator_filename = os.path.join(outputdir, 'dunelinearellipticcg2dsgrid.cc')
    module = pybindgen.Module('dunelinearellipticcg2dsgrid')
    module.add_include('"detailed_discretizations.hh"')
    #module.add_function('run', None, [])
    #namespace = module.add_cpp_namespace('Dune').add_cpp_namespace('Stuff').add_cpp_namespace('LA').add_cpp_namespace('Container')
    module.add_container('std::vector< std::string >', 'std::string', 'list')
    module.add_function('writeDescriptionFile', None, [pybindgen.param('std::string', 'filename'),
                                                       pybindgen.param('std::string', '_id')])
    MyVector = module.add_class('MyVector')
    MyVector.add_constructor([pybindgen.param('double', 's')])

    Operator = module.add_class('Operator')
    Operator.add_constructor([])
    Operator.add_method('apply',
                        pybindgen.retval('MyVector *', caller_owns_return=True),
                        [pybindgen.param('MyVector *', 'vector', transfer_ownership=False)])

    LinearEllipticExampleCG = module.add_class('LinearEllipticExampleCG')
    LinearEllipticExampleCG.add_constructor([pybindgen.param('std::vector< std::string >', 'arg')])
    LinearEllipticExampleCG.add_method('solve',
                                       pybindgen.retval('MyVector *', caller_owns_return=True),
                                       [])
    #EigenManager.add_method('createDenseVector',
                            #pybindgen.retval('std::string'),
                            #[pybindgen.param('int', 'size')])

    with open(generator_filename, 'wb') as output:
        module.generate(FileCodeSink(output))

    ## check if generate() was successfull
    #print os.stat(generator_fn).st_size

if __name__ == '__main__':
    generate_my_module(sys.argv[1], sys.argv[2], sys.argv[3:])
