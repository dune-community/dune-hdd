#ifndef DUNE_HDD_EXAMPLES_LINEARELLIPTIC_CG_MPI_HH
#define DUNE_HDD_EXAMPLES_LINEARELLIPTIC_CG_MPI_HH

#include <memory>
#include <dune/grid/spgrid.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/hdd/linearelliptic/discreteproblem.hh>
#include <dune/hdd/linearelliptic/discretizations/cg.hh>

class LinearellipticExampleCG
{
public:
  typedef Dune::SPGrid<double,2> GridType;
  typedef Dune::HDD::LinearElliptic::DiscreteProblem< GridType > DiscreteProblemType;
  typedef typename DiscreteProblemType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = DiscreteProblemType::dimRange;
  typedef Dune::HDD::LinearElliptic::Discretizations::CG
      < GridType, Dune::Stuff::Grid::ChooseLayer::leaf, RangeFieldType, dimRange, 1,
        Dune::GDT::ChooseSpaceBackend::pdelab, Dune::Stuff::LA::ChooseBackend::istl_sparse > DiscretizationType;

  static std::string static_id()
  {
    return "linearelliptic.cg";
  }

  static void write_config_file(const std::string filename = static_id() + ".cfg")
  {
    DiscreteProblemType::write_config(filename, static_id());
  }

  void initialize(const std::vector< std::string > arguments)
  {
    using namespace Dune;
    if (!discrete_problem_) {
      discrete_problem_ = Stuff::Common::make_unique< DiscreteProblemType >(static_id(), arguments);
      const bool debug_logging = discrete_problem_->debug_logging();
      auto& info = DSC_LOG_INFO;
      auto& debug = DSC_LOG_DEBUG;
      discretization_ = Stuff::Common::make_unique< DiscretizationType >(discrete_problem_->grid_provider(),
                                                                         discrete_problem_->boundary_info(),
                                                                         discrete_problem_->problem());
      info << "initializing discretization";
      if (debug_logging)
        info << ":" << std::endl;
      else
        info << "... " << std::flush;
      Dune::Timer timer;
      discretization_->init(debug, "  ");
      if (!debug_logging)
        info << "done (took " << timer.elapsed() << "s)" << std::endl;
    } else
      DSC_LOG_INFO << "initialize has already been called" << std::endl;
  } // ... initialize(...)

  const DiscreteProblemType& discrete_problem() const
  {
    return *discrete_problem_;
  }

  const DiscretizationType& discretization() const
  {
    return *discretization_;
  }

  DiscretizationType* discretization_and_return_ptr() const
  {
    return new DiscretizationType(*discretization_);
  }

private:
  std::unique_ptr< DiscreteProblemType > discrete_problem_;
  std::unique_ptr< DiscretizationType > discretization_;
}; // class LinearellipticExampleCG


#endif // CG_MPI_HH
