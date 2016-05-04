# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_python_bindings(linearellipticexamplegenericmultiscale generic_multiscale_bindings_generator.py generic_multiscale_bindings_generator.hh "")

foreach(target
        linearellipticexamplegenericmultiscale)
  add_dune_superlu_flags(${target})
  add_dune_alugrid_flags(${target})
  add_dune_tbb_flags(${target})
endforeach(target)

foreach(file
        "generic_multiscale.py")
  if (NOT EXISTS ${file})
    execute_process(COMMAND ln -s "${CMAKE_CURRENT_SOURCE_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${file}")
  endif (NOT EXISTS ${file})
endforeach(file)
