# This file is part of the dune-hdd project:
#   https://github.com/pymor/dune-hdd
# Copyright Holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

foreach(file
        "ORS2016__3_3.py"
        "morepas3.py"
        "morepas3__estimate.py"
        "morepas3__prepare.py"
        "morepas3__reduce.py")
  if (NOT EXISTS ${file})
    execute_process(COMMAND ln -s "${CMAKE_CURRENT_SOURCE_DIR}/${file}" "${CMAKE_CURRENT_BINARY_DIR}/${file}")
  endif (NOT EXISTS ${file})
endforeach(file)
