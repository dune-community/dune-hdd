# This module has been archived, no further development is expected.

    # This file is part of the dune-hdd project:
    #   http://users.dune-project.org/projects/dune-hdd
    # Copyright holders: Felix Schindler
    # License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

    dune-hdd is a DUNE (http://www.dune-project.org) module which provides high
    dimensional discretizations for grid-based discretizations of pdes. It makes
    use of the discretization modules dune-gdt
    (http://users.dune-project.org/projects/dune-gdt/) and dune-pdelab
    (http://www.dune-project.org/pdelab/index.html). For each pde, this module
    provides analytical problems (that describe the pde) together with
    appropriate discretizations for the pde. dune-hhd is also intended to be used
    with pyMor (http://pymor.org) in the context of model reduction, hence its
    name.

    New users may best try out this module by using the git supermodule
    dune-hdd-demos (http://users.dune-project.org/projects/dune-hdd-demos).
    Experienced DUNE users may go ahead. As usual, you will have to call autogen,
    configure and make using dunecontrol
    (see http://www.dune-project.org/doc/installation-notes.html), working
    examples are located in 'examples/'...
