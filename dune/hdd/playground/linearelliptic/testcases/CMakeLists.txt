# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

add_library(linearelliptic_testcase_esv2007_2daluconform
            OBJECT EXCLUDE_FROM_ALL
            ESV2007_2daluconform.cc)
add_library(linearelliptic_testcase_esv2007_multiscale_2daluconform
            OBJECT EXCLUDE_FROM_ALL
            ESV2007_multiscale_2daluconform.cc)
