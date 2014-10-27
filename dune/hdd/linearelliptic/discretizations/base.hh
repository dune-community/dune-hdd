// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH

#include <map>
#include <algorithm>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/common/logging.hh>

#include <dune/pymor/operators/base.hh>
#include <dune/pymor/operators/affine.hh>
#include <dune/pymor/functionals/affine.hh>

#include <dune/gdt/discretefunction/default.hh>

#include "interfaces.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward
template< class ImpTraits >
class ContainerBasedDefault;


namespace internal {


template< class MatrixImp, class VectorImp >
class ContainerBasedDefaultTraits
{
public:
  typedef MatrixImp MatrixType;
  typedef VectorImp VectorType;
  typedef Pymor::Operators::LinearAffinelyDecomposedContainerBased< MatrixType, VectorType > OperatorType;
  typedef OperatorType ProductType;
  typedef Pymor::Functionals::LinearAffinelyDecomposedVectorBased< VectorType > FunctionalType;
}; // class ContainerBasedDefaultTraits


} // namespace internal


template< class ImpTraits >
class CachedDefault
  : public DiscretizationInterface< ImpTraits >
{
  typedef DiscretizationInterface< ImpTraits >  BaseType;
  typedef CachedDefault< ImpTraits >            ThisType;
public:
  using typename BaseType::GridViewType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::ProblemType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::VectorType;

private:
  typedef Stuff::Grid::BoundaryInfoProvider< typename GridViewType::Intersection > BoundaryInfoProvider;

public:
  static std::string static_id() { return "hdd.linearelliptic.discretizations.cached"; }

  CachedDefault(const std::shared_ptr< const TestSpaceType >& test_spc,
                const std::shared_ptr< const AnsatzSpaceType > ansatz_spc,
                const Stuff::Common::Configuration& bnd_inf_cfg,
                const ProblemType& prb)
    : BaseType(prb)
    , test_space_(test_spc)
    , ansatz_space_(ansatz_spc)
    , boundary_info_cfg_(bnd_inf_cfg)
    , boundary_info_(BoundaryInfoProvider::create(boundary_info_cfg_.get< std::string >("type"), boundary_info_cfg_))
    , problem_(prb)
  {}

  CachedDefault(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  const std::shared_ptr< const TestSpaceType >& test_space() const
  {
    return test_space_;
  }

  const std::shared_ptr< const AnsatzSpaceType >& ansatz_space() const
  {
    return test_space_;
  }

  const GridViewType& grid_view() const
  {
    return test_space_->grid_view();
  }

  const Stuff::Common::Configuration& boundary_info_cfg() const
  {
    return boundary_info_cfg_;
  }

  const BoundaryInfoType& boundary_info() const
  {
    return *boundary_info_;
  }

  const ProblemType& problem() const
  {
    return problem_;
  }

  VectorType create_vector() const
  {
    return VectorType(ansatz_space_->mapper().size());
  }

  void visualize(const VectorType& vector,
                 const std::string filename,
                 const std::string name,
                 Pymor::Parameter mu = Pymor::Parameter()) const
  {
    VectorType tmp = vector.copy();
    const auto vectors = this->available_vectors();
    const auto result = std::find(vectors.begin(), vectors.end(), "dirichlet");
    if (result != vectors.end()) {
      const auto dirichlet_vector = this->get_vector("dirichlet");
      if (dirichlet_vector.parametric()) {
        const Pymor::Parameter mu_dirichlet = this->map_parameter(mu, "dirichlet");
        if (mu_dirichlet.type() != dirichlet_vector.parameter_type())
          DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                     mu_dirichlet.type() << " vs. " << dirichlet_vector.parameter_type());
        tmp = dirichlet_vector.freeze_parameter(mu);
      } else
        tmp = *(dirichlet_vector.affine_part());
      tmp += vector;
    }
    const GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType >function(*(ansatz_space()), tmp, name);
    function.visualize(filename);
  } // ... visualize(...)

  using BaseType::solve;

  void solve(const DSC::Configuration options, VectorType& vector, const Pymor::Parameter mu = Pymor::Parameter()) const
  {
    auto logger = DSC::TimedLogger().get(static_id());
    const auto options_in_cache = cache_.find(options);
    bool exists = (options_in_cache != cache_.end());
    typename std::map< Pymor::Parameter, std::shared_ptr< VectorType > >::const_iterator options_and_mu_in_cache;
    if (exists) {
      options_and_mu_in_cache = options_in_cache->second.find(mu);
      exists = (options_and_mu_in_cache != options_in_cache->second.end());
    }
    if (!exists) {
      logger.info() << "solving";
      if (!mu.empty())
        logger.info() << " for mu = " << mu;
      logger.info() << "... " << std::endl;
      uncached_solve(options, vector, mu);
      cache_[options][mu] = Stuff::Common::make_unique< VectorType >(vector.copy());
    } else {
      logger.info() << "retrieving solution ";
      if (!mu.empty())
        logger.info() << "for mu = " << mu << " ";
      logger.info() << "from cache... " << std::endl;
      const auto& result = *(options_and_mu_in_cache->second);
      vector = result;
    }
  } // ... solve(...)

  void uncached_solve(const DSC::Configuration options, VectorType& vector, const Pymor::Parameter mu) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().uncached_solve(options, vector, mu));
  }

protected:
  const std::shared_ptr< const TestSpaceType > test_space_;
  const std::shared_ptr< const AnsatzSpaceType > ansatz_space_;
  const Stuff::Common::Configuration boundary_info_cfg_;
  const std::shared_ptr< const BoundaryInfoType > boundary_info_;
  const ProblemType& problem_;

  mutable std::map< DSC::Configuration, std::map< Pymor::Parameter, std::shared_ptr< VectorType > > > cache_;
}; // class CachedDefault


template< class ImpTraits >
class ContainerBasedDefault
  : public CachedDefault< ImpTraits >
{
  typedef CachedDefault< ImpTraits >         BaseType;
  typedef ContainerBasedDefault< ImpTraits > ThisType;
public:
  typedef ImpTraits Traits;
  typedef typename Traits::MatrixType     MatrixType;
  typedef typename Traits::OperatorType   OperatorType;
  typedef typename Traits::ProductType    ProductType;
  typedef typename Traits::FunctionalType FunctionalType;

  using typename BaseType::GridViewType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::ProblemType;
  using typename BaseType::TestSpaceType;
  using typename BaseType::AnsatzSpaceType;
  using typename BaseType::VectorType;

protected:
  typedef Pymor::LA::AffinelyDecomposedContainer< MatrixType > AffinelyDecomposedMatrixType;
  typedef Pymor::LA::AffinelyDecomposedContainer< VectorType > AffinelyDecomposedVectorType;

public:
  static std::string static_id() { return "hdd.linearelliptic.discretizations.containerbased"; }

  ContainerBasedDefault(const std::shared_ptr< const TestSpaceType >& test_spc,
                        const std::shared_ptr< const AnsatzSpaceType >& ansatz_spc,
                        const Stuff::Common::Configuration& bnd_inf_cfg,
                        const ProblemType& prb)
    : BaseType(test_spc, ansatz_spc, bnd_inf_cfg, prb)
    , container_based_initialized_(false)
    , purely_neumann_(false)
    , matrix_(std::make_shared< AffinelyDecomposedMatrixType >())
    , rhs_(std::make_shared< AffinelyDecomposedVectorType >())
  {}

  ContainerBasedDefault(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  std::shared_ptr< AffinelyDecomposedMatrixType >& system_matrix()
  {
    return matrix_;
  }

  const std::shared_ptr< const AffinelyDecomposedMatrixType >& system_matrix() const
  {
    return matrix_;
  }

  std::shared_ptr< AffinelyDecomposedVectorType > rhs()
  {
    return rhs_;
  }

  const std::shared_ptr< const AffinelyDecomposedVectorType >& rhs() const
  {
    return rhs_;
  }

  OperatorType get_operator() const
  {
    assert_everything_is_ready();
    return OperatorType(*matrix_);
  }

  FunctionalType get_rhs() const
  {
    assert_everything_is_ready();
    return FunctionalType(*rhs_);
  }

  std::vector< std::string > available_products() const
  {
    if (products_.size() == 0)
      return std::vector< std::string >();
    std::vector< std::string > ret;
    for (const auto& pair : products_)
      ret.push_back(pair.first);
    return ret;
  } // ... available_products(...)

  ProductType get_product(const std::string id) const
  {
    if (products_.size() == 0)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Do not call get_product() if available_products() is empty!");
    const auto result = products_.find(id);
    if (result == products_.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, id);
    return ProductType(*(result->second));
  } // ... get_product(...)

  std::vector< std::string > available_vectors() const
  {
    if (vectors_.size() == 0)
      return std::vector< std::string >();
    std::vector< std::string > ret;
    for (const auto& pair : vectors_)
      ret.push_back(pair.first);
    return ret;
  } // ... available_vectors(...)

  AffinelyDecomposedVectorType get_vector(const std::string id) const
  {
    if (vectors_.size() == 0)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "Do not call get_vector() if available_vectors() is empty!");
    const auto result = vectors_.find(id);
    if (result == vectors_.end())
      DUNE_THROW(Stuff::Exceptions::wrong_input_given, id);
    return *(result->second);
  } // ... get_vector(...)

  /**
   * \brief solves for u_0
   */
  void uncached_solve(VectorType& vector, const Pymor::Parameter mu = Pymor::Parameter()) const
  {
    auto logger = DSC::TimedLogger().get(static_id());
    assert_everything_is_ready();
    if (mu.type() != this->parameter_type())
      DUNE_THROW(Pymor::Exceptions::wrong_parameter_type, mu.type() << " vs. " << this->parameter_type());
    const auto& rhs = *(this->rhs_);
    const auto& matrix = *(this->matrix_);
    if (purely_neumann_) {
      VectorType tmp_rhs = rhs.parametric() ? rhs.freeze_parameter(this->map_parameter(mu, "rhs"))
                                            : *(rhs.affine_part());
      MatrixType tmp_system_matrix = matrix.parametric() ? matrix.freeze_parameter(this->map_parameter(mu, "lhs"))
                                                         : *(matrix.affine_part());
      tmp_system_matrix.unit_row(0);
      tmp_rhs.set_entry(0, 0.0);
      Stuff::LA::Solver< MatrixType >(tmp_system_matrix).apply(tmp_rhs, vector);
      vector -= vector.mean();
    } else {
      // compute right hand side vector
      logger.debug() << "computing right hand side..." << std::endl;
      std::shared_ptr< const VectorType > rhs_vector;
      if (!rhs.parametric())
        rhs_vector = rhs.affine_part();
      else {
        const Pymor::Parameter mu_rhs = this->map_parameter(mu, "rhs");
        rhs_vector = std::make_shared< const VectorType >(rhs.freeze_parameter(mu_rhs));
      }
      logger.debug() << "computing system matrix..." << std::endl;
      const OperatorType lhsOperator(matrix);
      if (lhsOperator.parametric()) {
        const Pymor::Parameter mu_lhs = this->map_parameter(mu, "lhs");
        const auto frozenOperator = lhsOperator.freeze_parameter(mu_lhs);
        frozenOperator.apply_inverse(*rhs_vector, vector);
      } else {
        const auto nonparametricOperator = lhsOperator.affine_part();
        nonparametricOperator.apply_inverse(*rhs_vector, vector);
      }
    }
  } // ... uncached_solve(...)

protected:
  void assert_everything_is_ready() const
  {
    if (!container_based_initialized_)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "The implemented discretization has to fill 'matrix_' and 'rhs_' during init() and set "
                 << "container_based_initialized_ to true!\n"
                 << "The user has to call init() before calling any other method!");
  } // ... assert_everything_is_ready()

  bool container_based_initialized_;
  bool purely_neumann_;
  std::shared_ptr< AffinelyDecomposedMatrixType > matrix_;
  std::shared_ptr< AffinelyDecomposedVectorType > rhs_;
  mutable std::map< std::string, std::shared_ptr< AffinelyDecomposedMatrixType > > products_;
  mutable std::map< std::string, std::shared_ptr< AffinelyDecomposedVectorType > > vectors_;
}; // class ContainerBasedDefault


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH
