
SUBDIRS =

noinst_PROGRAMS = pdelab

pdelab_SOURCES = pdelab.cc

GRIDDIM=2
GRIDTYPE=SGRID

CUSTOM_CPPFLAGS = -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -Wall -O0 -DDEBUG -g -fdiagnostics-show-option -std=c++0x

pdelab_CPPFLAGS = $(AM_CPPFLAGS) \
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
pdelab_LDADD = \
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
pdelab_LDFLAGS = $(AM_LDFLAGS) \
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

include $(top_srcdir)/am/global-rules

