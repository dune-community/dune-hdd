# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-istl dune-fem dune-fem-tools dune-detailed-discretizations
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-fem], [dune-fem-tools], [dune-detailed-discretizations]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/detailed-solvers/Makefile
  dune/detailed-solvers/linearelliptic/Makefile
  dune/detailed-solvers/linearelliptic/model/Makefile
  dune/detailed-solvers/linearelliptic/model/affineparametric/Makefile
  dune/detailed-solvers/linearelliptic/solver/cg/Makefile
  dune/detailed-solvers/linearelliptic/solver/fv/Makefile
  dune/detailed-solvers/linearelliptic/solver/multiscale/Makefile
  dune/detailed-solvers/linearelliptic/solver/multiscale/semicg/Makefile
  m4/Makefile
  dune-detailed-solvers.pc
  examples/Makefile
  examples/linearelliptic/Makefile
  examples/linearelliptic/cg/Makefile
  examples/linearelliptic/fv/Makefile
  examples/linearelliptic/multiscale/Makefile
  examples/linearelliptic/multiscale/semicg/Makefile
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
