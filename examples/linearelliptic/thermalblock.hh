// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
#endif

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif


#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/provider/cube.hh>

#include <dune/hdd/linearelliptic/problems/thermalblock.hh>
#include <dune/hdd/linearelliptic/discretizations/cg.hh>

namespace internal {


template< class GridType >
class Initializer
{
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;

public:
  Initializer(const std::string& num_grid_elements,
              const DUNE_STUFF_SSIZE_T info_log_levels,
              const DUNE_STUFF_SSIZE_T debug_log_levels,
              const bool enable_warnings,
              const bool enable_colors,
              const std::string info_color,
              const std::string debug_color,
              const std::string warn_color)
  {
    try {
      int argc = 0;
      char** argv = new char* [0];
#if HAVE_DUNE_FEM
      Dune::Fem::MPIManager::initialize(argc, argv);
#else
      Dune::MPIHelper::instance(argc, argv);
#endif
    } catch (...) {}
    try {
      DSC::TimedLogger().create(info_log_levels,
                                debug_log_levels,
                                enable_warnings,
                                enable_colors,
                                info_color,
                                debug_color,
                                warn_color);
    } catch (Dune::Stuff::Exceptions::you_are_using_this_wrong&) {}
    DSC::TimedLogger().get("cg.thermalblock.example").info() << "creating grid and problem... " << std::endl;
    grid_provider_ = GridProviderType::create(GridProviderType::default_config().add(DSC::Configuration("num_elements",
                                                                                                        num_grid_elements),
                                                                                     "",
                                                                                     true));
#if HAVE_ALUGRID
    if (std::is_same< GridType, Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > >::value)
      grid_provider_->grid().globalRefine(1);
#endif // HAVE_ALUGRID
  } // Initializer(...)

protected:
  std::unique_ptr< GridProviderType > grid_provider_;
}; // class Initializer


} // namespace internal


template< class GridImp, Dune::GDT::ChooseSpaceBackend space_backend, Dune::Stuff::LA::ChooseBackend la_backend >
class CgExample
  : internal::Initializer< GridImp >
{
  typedef internal::Initializer< GridImp >               BaseType;
  typedef GridImp                                        GridType;
  typedef typename GridType::template Codim< 0 >::Entity EntityType;
  static const size_t                                    dimDomain = GridType::dimension;
  typedef typename GridType::ctype                       DomainFieldType;
  typedef double                                         RangeFieldType;
  static const size_t                                    dimRange = 1;
  typedef Dune::HDD::LinearElliptic::Problems::Thermalblock
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 > ProblemType;
public:
  typedef Dune::HDD::LinearElliptic::Discretizations::CG< GridType, Dune::Stuff::Grid::ChooseLayer::leaf,
                                                          RangeFieldType, dimRange, 1,
                                                          space_backend, la_backend > DiscretizationType;

  CgExample(const std::string& num_blocks        = "[2 2 2]",
            const std::string& num_grid_elements = "[32 32 32]",
            const DUNE_STUFF_SSIZE_T info_log_levels  = 0,
            const DUNE_STUFF_SSIZE_T debug_log_levels = -1,
            const bool enable_warnings = true,
            const bool enable_colors   = true,
            const std::string info_color  = DSC::TimedLogging::default_info_color(),
            const std::string debug_color = DSC::TimedLogging::default_debug_color(),
            const std::string warn_color  = DSC::TimedLogging::default_warning_color())
    : BaseType(num_grid_elements,
               info_log_levels,
               debug_log_levels,
               enable_warnings,
               enable_colors,
               info_color,
               debug_color,
               warn_color)
    , problem_(ProblemType::create(ProblemType::default_config().add(DSC::Configuration("diffusion_factor.num_elements",
                                                                                        num_blocks),
                                                                     "",
                                                                     true)))
    , discretization_(*grid_provider_, Dune::Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config(), *problem_)
  {
    auto logger = DSC::TimedLogger().get("cg.thermalblock.example");
    logger.info() << "initializing discretization... " << std::flush;
    discretization_.init();
    logger.info() << "done (grid has " << discretization_.grid_view().indexSet().size(0)
                  << " elements, discretization has " << discretization_.ansatz_space().mapper().size() << " DoFs)"
                  << std::endl;
  } // ... CgExample(...)

  DiscretizationType& discretization()
  {
    return discretization_;
  }

private:
  using BaseType::grid_provider_;
  std::unique_ptr< ProblemType > problem_;
  DiscretizationType discretization_;
}; // class CgExample


/**
    Orthonormalize a |VectorArray| using the stabilized Gram-Schmidt algorithm.

    Parameters
    ----------
    A
        The |VectorArray| which is to be orthonormalized.
    product
        The scalar product w.r.t. which to orthonormalize, given as a linear
        |Operator|. If `None` the Euclidean product is used.
    atol
        Vectors of norm smaller than `atol` are removed from the array.
    rtol
        Relative tolerance used to detect linear dependent vectors
        (which are then removed from the array).
    offset
        Assume that the first `offset` vectors are already orthogonal and start the
        algorithm at the `offset + 1`-th vector.
    find_duplicates
        If `True`, eliminate duplicate vectors before the main loop.
    reiterate
        If `True`, orthonormalize again if the norm of the orthogonalized vector is
        much smaller than the norm of the original vector.
    reiteration_threshold
        If `reiterate` is `True`, re-orthonormalize if the ratio between the norms of
        the orthogonalized vector and the original vector is smaller than this value.
    check
        If `True`, check if the resulting VectorArray is really orthonormal.
    check_tol
        Tolerance for the check.
    copy
        If `True`, create a copy of `A` instead of modifying `A` itself.


    Returns
    -------
    The orthonormalized |VectorArray|.
  **/
template< class VectorType >
std::vector< VectorType > gram_schmidt(const std::vector< VectorType >& A,
                                       const bool reiterate = true,
                                       const bool check = true,
                                       const double atol = 1e-13,
                                       const double rtol = 1e-13,
                                       const DUNE_STUFF_SSIZE_T offset = 0,
                                       const double reiteration_threshold = 1e-1,
                                       const double check_tol = 1e-3)
{
  if (offset < 0)
    DUNE_THROW(Dune::Stuff::Exceptions::index_out_of_range, offset);
  auto logger = DSC::TimedLogger().get("gram_schmidt");

  auto ret = A;

  // main loop
  std::set< size_t > remove;
  for (size_t i = offset; i < ret.size(); ++i) {
    // first calculate norm
    const double initial_norm = ret[i].l2_norm();

    if (initial_norm < atol) {
      logger.info() << "Removing vector " << i << " of norm " << initial_norm << std::endl;
      remove.insert(i);
    }

    if (i == 0)
      ret[0].scal(1./initial_norm);
    else {
      bool first_iteration = true;
      double norm = initial_norm;
      double old_norm = initial_norm;
      // If reiterate is true, reiterate as long as the norm of the vector changes
      // strongly during orthonormalization (due to Andreas Buhr).
      while (first_iteration || reiterate && norm/old_norm < reiteration_threshold) {

        if (first_iteration)
          first_iteration = false;
        else
          logger.info() << "Orthonormalizing vector " << i << " again" << std::endl;

        // orthogonalize to all vectors left
        for (size_t j = 0; j < i; ++j) {
          if (remove.find(j) != remove.end())
            continue;
          const double p = ret[i].dot(ret[j]);
          ret[i].axpy(-p, ret[j]);
        }

        // calculate new norm
        old_norm = norm;
        norm = ret[i].l2_norm();

        // remove vector if it got too small:
        if (norm/initial_norm < rtol) {
          logger.info() << "Removing linear dependent vector " << i << std::endl;
          remove.insert(i);
          break;
        }
      }
      ret[i].scal(1./norm);
    }
  }

  for (const size_t& i : remove)
    ret.erase(ret.begin() + i);

  if (check) {
    for (size_t i = offset; i < ret.size(); ++i)
      for (size_t j = i; j < ret.size(); ++j)
        if (std::abs(ret[i].dot(ret[j])) - (i == j ? 1. : 0.) > check_tol)
          DUNE_THROW(Dune::MathError,
                     "result not orthogonal: \n"
                     << "  A[i]*A[j] = " << ret[i].dot(ret[j]) << "\n"
                     << "  i = " << i << "\n"
                     << "  j = " << j);
  }
  return ret;
} // ... gram_schmidt(...)

template< class VectorType, class MatrixType >
std::vector< VectorType > gram_schmidt(const std::vector< VectorType >& A,
                                       const MatrixType& product,
                                       const bool reiterate = true,
                                       const bool check = true,
                                       const double atol = 1e-13,
                                       const double rtol = 1e-13,
                                       const DUNE_STUFF_SSIZE_T offset = 0,
                                       const double reiteration_threshold = 1e-1,
                                       const double check_tol = 1e-3)
{
  if (offset < 0)
    DUNE_THROW(Dune::Stuff::Exceptions::index_out_of_range, offset);
  auto logger = DSC::TimedLogger().get("gram_schmidt");

  auto ret = A;
  VectorType tmp(product.rows());

  // main loop
  std::set< size_t > remove;
  for (size_t i = offset; i < ret.size(); ++i) {
    // first calculate norm
    product.mv(ret[i], tmp);
    const double initial_norm = std::sqrt(tmp*ret[i]);

    if (initial_norm < atol) {
      logger.info() << "Removing vector " << i << " of norm " << initial_norm << std::endl;
      remove.insert(i);
    }

    if (i == 0)
      ret[0].scal(1./initial_norm);
    else {
      bool first_iteration = true;
      double norm = initial_norm;
      double old_norm = initial_norm;
      // If reiterate is true, reiterate as long as the norm of the vector changes
      // strongly during orthonormalization (due to Andreas Buhr).
      while (first_iteration || reiterate && norm/old_norm < reiteration_threshold) {

        if (first_iteration)
          first_iteration = false;
        else
          logger.info() << "Orthonormalizing vector " << i << " again" << std::endl;

        // orthogonalize to all vectors left
        for (size_t j = 0; j < i; ++j) {
          if (remove.find(j) != remove.end())
            continue;
          product.mv(ret[i], tmp);
          const double p = tmp*ret[j];
          ret[i].axpy(-p, ret[j]);
        }

        // calculate new norm
        old_norm = norm;
        product.mv(ret[i], tmp);
        norm = std::sqrt(tmp*ret[i]);

        // remove vector if it got too small:
        if (norm/initial_norm < rtol) {
          logger.info() << "Removing linear dependent vector " << i << std::endl;
          remove.insert(i);
          break;
        }
      }
      ret[i].scal(1./norm);
    }
  }

  for (const size_t& i : remove)
    ret.erase(ret.begin() + i);

  if (check) {
    for (size_t i = offset; i < ret.size(); ++i) {
      product.mv(ret[i], tmp);
      for (size_t j = i; j < ret.size(); ++j) {
        if (tmp*ret[j] - (i == j ? 1. : 0.) > check_tol)
          DUNE_THROW(Dune::MathError,
                     "result not orthogonal: \n"
                     << "  product.apply2(A[i], A[j]) = " << tmp*ret[j] << "\n"
                     << "  i = " << i << "\n"
                     << "  j = " << j);
      }
    }
  }
  return ret;
} // ... gram_schmidt(...)

#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_THERMALBLOCK_HH
