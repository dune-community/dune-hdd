#! /usr/bin/env python

# This file is part of the dune-pymor project:
#   https://github.com/pyMor/dune-pymor
# Copyright Holders: Felix Albrecht, Stephan Rave
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import sys, os
import pybindgen
from dune.pymor.core import create_module as create_dune_module
from dune.pymor.la.container import add_VectorInterface as add_dune_Vector


if __name__ == '__main__':
    # prepare
    assert(len(sys.argv) > 1)
    inputdir, outputdir, includedirs = sys.argv[1], sys.argv[2], sys.argv[3:]
    includedirs = includedirs[0].split(';')
    global_name = 'lacontainerexample'
    generator_filename = os.path.join(outputdir, global_name + '_bindings_generator.cc')
    # add dune-pymor constructs
    module, exceptions = create_dune_module(global_name, 'container.hh', add_all_of_dune_pymor=False)
    module, _ = add_dune_Vector(module, exceptions, 'Dune::Pymor::LA::DuneDynamicVector')
    module, _ = add_dune_Vector(module, exceptions, 'Dune::Pymor::LA::EigenDenseVector')
    # add user code stuff
    module.add_function('createDuneDynamicVector',
                        pybindgen.retval('Dune::Pymor::LA::DuneDynamicVector *', caller_owns_return=True),
                        [pybindgen.param('int', 'ss')])
    module.add_function('createEigenDenseVector',
                        pybindgen.retval('Dune::Pymor::LA::EigenDenseVector *', caller_owns_return=True),
                        [pybindgen.param('int', 'ss')])
    with open(generator_filename, 'wb') as output:
        module.generate(pybindgen.FileCodeSink(output))
