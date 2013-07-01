# This file is part of the dune-hdd project:
#   http://users.dune-project.org/projects/dune-hdd
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

SUBDIRS =

noinst_PROGRAMS = examples_linearelliptic_cg_gdt

examples_linearelliptic_cg_gdt_SOURCES = gdt.cc

GRIDDIM=2
GRIDTYPE=ALUGRID_SIMPLEX

CUSTOM_CPPFLAGS = -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -Wall -O0 -DDEBUG -g -fdiagnostics-show-option -std=c++0x

examples_linearelliptic_cg_gdt_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS) \
  $(EIGEN_CFLAGS) \
  $(BOOST_CPPFLAGS) \
  $(CUSTOM_CPPFLAGS)
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
examples_linearelliptic_cg_gdt_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD) \
  $(BOOST_LDADD) \
  -lboost_system \
  -lboost_filesystem
examples_linearelliptic_cg_gdt_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS) \
  $(BOOST_LDFLAGS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-helper-tools=$(DUNE_FEM_TOOLS_ROOT) --with-dune-gdt=$(DUNE_GDT_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

