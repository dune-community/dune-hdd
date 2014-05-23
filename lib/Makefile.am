# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

dunemodulelibdir=$(libdir)

#the dune-hdd library
lib_LTLIBRARIES = libdunehdd.la

# this is just a renaming of libhdd to prevent name clashes
libdunehdd_la_SOURCES =
# This forces automake to use the C++ linker
# (see the automake manual, section "Libtool Convenience Libraries")
nodist_EXTRA_libdunehdd_la_SOURCES = dummy.cc
sourcescheck_DUMMY = dummy.cc
libdunehdd_la_LIBADD = ../dune/hdd/libhdd.la

include $(top_srcdir)/am/global-rules
