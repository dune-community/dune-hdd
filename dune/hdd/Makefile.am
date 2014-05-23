This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

hddincludedir = $(includedir)/dune/hdd/
hddinclude_HEADERS =

SUBDIRDS = .

noinst_LTLIBRARIES = libhdd.la

libhdd_la_SOURCES = \
  linearelliptic/problems.cc \
  linearelliptic/discreteproblem.cc \
  linearelliptic/discretizations/cg_alugrid_conform_2d_fem.cc \
  linearelliptic/discretizations/cg_alugrid_conform_2d_pdelab.cc \
  linearelliptic/discretizations/swipdg_alugrid_conform_2d_femlocalfunctions.cc

libhdd_la_LIBADD = $(DUNE_LIBS) $(ALUGRID_LIBS)

libhdd_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(ALUGRID_CPPFLAGS)

include $(top_srcdir)/am/global-rules
