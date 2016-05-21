#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_MULTISCALE_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_MULTISCALE_HH

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
# include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/container-interface.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/grid/multiscale/provider.hh>

#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/spaces/interface.hh>

#include <dune/hdd/linearelliptic/discretizations/block-swipdg.hh>
#include <dune/hdd/linearelliptic/problems.hh>


template< class GridImp, Dune::GDT::ChooseSpaceBackend space_backend, Dune::Stuff::LA::ChooseBackend la_backend >
class GenericLinearellipticMultiscaleExample
{
  typedef GridImp GridType;
  typedef typename GridType::template Codim< 0 >::Entity E;
  typedef typename GridType::ctype D;
  static const size_t d = GridType::dimension;
  typedef double R;
  static const size_t r = 1;
public:
  typedef Dune::HDD::LinearElliptic::Discretizations::BlockSWIPDG< GridType, R, 1, 1, la_backend > DiscretizationType;
  typedef typename DiscretizationType::VectorType                                     VectorType;
private:
  typedef typename DiscretizationType::MatrixType                                        MatrixType;
  typedef Dune::grid::Multiscale::MsGridProviders< GridType >                            GridProvider;
  typedef Dune::Stuff::Grid::BoundaryInfoProvider< typename GridType::LeafIntersection > BoundaryProvider;
  typedef Dune::HDD::LinearElliptic::ProblemsProvider< E, D, d, R, r >                   ProblemProvider;
  typedef Dune::Stuff::LA::Solver< MatrixType >                                          SolverProvider;
public:
  static Dune::Stuff::Common::Configuration logger_options();

  static std::vector< std::string > grid_options();

  static Dune::Stuff::Common::Configuration grid_options(const std::string& type);

  static std::vector< std::string > boundary_options();

  static Dune::Stuff::Common::Configuration boundary_options(const std::string& type);

  static std::vector< std::string > problem_options();

  static Dune::Stuff::Common::Configuration problem_options(const std::string& type);

  static std::vector< std::string > solver_options();

  static Dune::Stuff::Common::Configuration solver_options(const std::string& type);

  GenericLinearellipticMultiscaleExample(const Dune::Stuff::Common::Configuration& logger_cfg = Dune::Stuff::Common::Configuration(),
                                         const Dune::Stuff::Common::Configuration& grid_cfg = Dune::Stuff::Common::Configuration(),
                                         const Dune::Stuff::Common::Configuration& boundary_cfg = Dune::Stuff::Common::Configuration(),
                                         const Dune::Stuff::Common::Configuration& problem_cfg = Dune::Stuff::Common::Configuration());

  DiscretizationType& discretization();

  void visualize(const std::string& filename_prefix) const;

  VectorType project(const std::string& expression) const;

  VectorType prolong(const DiscretizationType& source_disc, const VectorType& source_vec) const;

  VectorType oswald_interpolate(const VectorType& vector) const;

  double alpha(const Dune::Pymor::Parameter& mu_1, const Dune::Pymor::Parameter& mu_2) const;

  double gamma(const Dune::Pymor::Parameter& mu_1, const Dune::Pymor::Parameter& mu_2) const;

  double min_diffusion_ev(const Dune::Pymor::Parameter& mu) const;

  double max_diffusion_ev(const Dune::Pymor::Parameter& mu) const;

  double elliptic_reconstruction_estimate(const VectorType& p_h,
                                          const Dune::Pymor::Parameter& mu_min,
                                          const Dune::Pymor::Parameter& mu_max,
                                          const Dune::Pymor::Parameter& mu_hat,
                                          const Dune::Pymor::Parameter& mu_bar,
                                          const Dune::Pymor::Parameter& mu) const;

private:
  const Dune::Stuff::Common::Configuration boundary_cfg_;
  std::unique_ptr< Dune::grid::Multiscale::ProviderInterface< GridType > > grid_;
  std::unique_ptr< Dune::HDD::LinearElliptic::ProblemInterface< E, D, d, R, r > > problem_;
  std::unique_ptr< DiscretizationType > discretization_;

}; // class GenericLinearellipticMultiscaleExample


#if HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL


extern template class GenericLinearellipticMultiscaleExample< Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming >,
                                                              Dune::GDT::ChooseSpaceBackend::fem,
                                                              Dune::Stuff::LA::ChooseBackend::istl_sparse >;


#endif // HAVE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL


#endif // DUNE_HDD_EXAMPLES_LINEARELLIPTIC_GENERIC_MULTISCALE_HH
