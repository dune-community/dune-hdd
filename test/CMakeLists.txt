# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

enable_testing()

DEPENDENCYCHECK(${hdd})

BEGIN_TESTCASES()
END_TESTCASES()

add_dune_alugrid_flags(test_OS2014)
