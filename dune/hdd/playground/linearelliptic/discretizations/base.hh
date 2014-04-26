// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH

#include <map>
#include <algorithm>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/exceptions.hh>

#include <dune/pymor/operators/base.hh>
#include <dune/pymor/operators/affine.hh>
#include <dune/pymor/functionals/affine.hh>

#include <dune/gdt/discretefunction/default.hh>

#include "../../../linearelliptic/discretizations/interfaces.hh"

namespace Dune {
namespace HDD {
namespace LinearElliptic {
namespace Discretizations {


// forward
template< class ImpTraits >
class ContainerBasedDefault;


namespace internal {


template< class ImpTraits >
class ContainerBasedDefaultTraits
  : public ImpTraits
{
public:
  typedef ContainerBasedDefault< ImpTraits > derived_type;
  using typename ImpTraits::MatrixType;
  using typename ImpTraits::VectorType;
  typedef Pymor::Operators::LinearAffinelyDecomposedContainerBased< MatrixType, VectorType > OperatorType;
  typedef OperatorType ProductType;
  typedef Pymor::Functionals::LinearAffinelyDecomposedVectorBased< VectorType > FunctionalType;
}; // class ContainerBasedDefaultTraits


} // namespace internal


template< class ImpTraits >
class CachedDefault
  : public DiscretizationInterface< ImpTraits >
{
  typedef DiscretizationInterface< ImpTraits > BaseType;
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
  CachedDefault(const std::shared_ptr< const TestSpaceType >& test_spc,
                const std::shared_ptr< const AnsatzSpaceType > ansatz_spc,
                const Stuff::Common::ConfigTree& bnd_inf_cfg,
                const ProblemType& prb)
    : BaseType(prb)
    , test_space_(test_spc)
    , ansatz_space_(ansatz_spc)
    , boundary_info_cfg_(bnd_inf_cfg)
    , boundary_info_(BoundaryInfoProvider::create(boundary_info_cfg_.get< std::string >("type"), boundary_info_cfg_))
    , problem_(prb)
  {}

  CachedDefault(const CachedDefault& other) = delete;

  CachedDefault& operator=(const CachedDefault& other) = delete;

  const std::shared_ptr< const TestSpaceType >& test_space() const
  {
    return test_space_;
  }

  const std::shared_ptr< const AnsatzSpaceType >& ansatz_space() const
  {
    return test_space_;
  }

  const std::shared_ptr< const GridViewType >& grid_view() const
  {
    return test_space_->grid_view();
  }

  const Stuff::Common::ConfigTree& boundary_info_cfg() const
  {
    return boundary_info_cfg_;
  }

  const BoundaryInfoType& boundary_info() const
  {
    return *boundary_info_;
  }

  const ProblemType& problem() const
  {
    return *problem_;
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
          DUNE_THROW_COLORFULLY(Pymor::Exceptions::wrong_parameter_type,
                                mu_dirichlet.type() << " vs. " << dirichlet_vector.parameter_type());
        tmp = dirichlet_vector.freeze_parameter(mu);
      } else
        tmp = *(dirichlet_vector.affine_part());
      tmp += vector;
    }
    const GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType >function(*(ansatz_space()), tmp, name);
    function.visualize(filename);
  } // ... visualize(...)

  void solve(VectorType& vector, const Pymor::Parameter mu = Pymor::Parameter()) const
  {
    const auto search_result = cache_.find(mu);
    if (search_result == cache_.end()) {
      uncached_solve(vector, mu);
      cache_.insert(std::make_pair(mu, Stuff::Common::make_unique< VectorType >(vector.copy())));
    } else {
      const auto& result = *(search_result->second);
      vector = result;
    }
  } // ... solve(...)

  void uncached_solve(VectorType& vector, const Pymor::Parameter mu) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).uncached_solve(vector, mu));
  }

protected:
  const std::shared_ptr< const TestSpaceType > test_space_;
  const std::shared_ptr< const AnsatzSpaceType > ansatz_space_;
  const Stuff::Common::ConfigTree& boundary_info_cfg_;
  const std::unique_ptr< const BoundaryInfoType > boundary_info_;
  const ProblemType& problem_;

  mutable std::map< Pymor::Parameter, std::unique_ptr< VectorType > > cache_;
}; // class CachedDefault


template< class ImpTraits >
class ContainerBasedDefault
  : public CachedDefault< internal::ContainerBasedDefaultTraits< ImpTraits > >
{
  typedef CachedDefault< internal::ContainerBasedDefaultTraits< ImpTraits > > BaseType;
public:
  typedef internal::ContainerBasedDefaultTraits< ImpTraits > Traits;
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
  ContainerBasedDefault(const std::shared_ptr< const TestSpaceType >& test_spc,
                        const std::shared_ptr< const AnsatzSpaceType > ansatz_spc,
                        const Stuff::Common::ConfigTree& bnd_inf_cfg,
                        const ProblemType& prb)
    : BaseType(test_spc, ansatz_spc, bnd_inf_cfg, prb)
    , matrix_(nullptr)
    , rhs_(nullptr)
  {}

  OperatorType get_operator() const
  {
    assert_everything_is_ready();
    return OperatorType(matrix_);
  }

  FunctionalType get_rhs() const
  {
    assert_everything_is_ready();
    return FunctionalType(rhs_);
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
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::you_are_using_this_wrongly,
                            "Do not call get_product() if available_products() is empty!");
    const auto result = products_.find(id);
    if (result == products_.end())
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given, id);
    return ProductType(result->second);
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
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::you_are_using_this_wrongly,
                            "Do not call get_vector() if available_vectors() is empty!");
    const auto result = vectors_.find(id);
    if (result == vectors_.end())
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::wrong_input_given, id);
    return *(result->second);
  } // ... get_vector(...)

  /**
   * \brief solves for u_0
   */
  void uncached_solve(VectorType& vector, const Pymor::Parameter mu = Pymor::Parameter()) const
  {
    if (mu.type() != this->parameter_type())
      DUNE_THROW_COLORFULLY(Pymor::Exceptions::wrong_parameter_type,
                            mu.type() << " vs. " << this->parameter_type());
//    Dune::Timer timer;
    // compute right hand side vector
    const auto& rhs = *(this->rhs_);
    std::shared_ptr< const VectorType > rhs_vector;
    if (!rhs.parametric())
      rhs_vector = rhs.affine_part();
    else {
//      out << prefix << "computing rhs... " << std::flush;
//      timer.reset();
      const Pymor::Parameter muRhs = this->map_parameter(mu, "rhs");
      rhs_vector = std::make_shared< const VectorType >(rhs.freeze_parameter(muRhs));
//      out << "done (took " << timer.elapsed() << "s)" << std::endl;
    }
    const OperatorType lhsOperator(*(this->matrix_));
    if (lhsOperator.parametric()) {
//      out << prefix << "computing lhs... " << std::flush;
//      timer.reset();
      const Pymor::Parameter muLhs = Pymor::Parametric::map_parameter(mu, "lhs");
      const auto frozenOperator = lhsOperator.freeze_parameter(muLhs);
//      out << "done (took " << timer.elapsed() << "s)" << std::endl;
//      const std::string option = frozenOperator.invert_options()[0];
//      out << prefix << "solving with '" << option << "' option... " << std::flush;
//      timer.reset();
      frozenOperator.apply_inverse(*rhs_vector, vector/*, option*/);
//      out << "done (took " << timer.elapsed() << "s)" << std::endl;
    } else {
      const auto nonparametricOperator = lhsOperator.affine_part();
      const std::string option = nonparametricOperator.invert_options()[0];
//      out << prefix << "solving with '" << option << "' option... " << std::flush;
//      timer.reset();
      nonparametricOperator.apply_inverse(*rhs_vector, vector, option);
//      out << "done (took " << timer.elapsed() << "s)" << std::endl;
    }
  } // ... uncached_solve(...)

protected:
  bool assert_everything_is_ready() const
  {
    if (!matrix_ || !rhs_)
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::internal_error,
                            "The discretization has to fill 'matrix_' and 'rhs_' during init()!");
  } // ... assert_everything_is_ready()

  mutable std::shared_ptr< AffinelyDecomposedMatrixType > matrix_;
  mutable std::shared_ptr< AffinelyDecomposedVectorType > rhs_;
  mutable std::map< std::string, std::shared_ptr< AffinelyDecomposedMatrixType > > products_;
  mutable std::map< std::string, std::shared_ptr< AffinelyDecomposedVectorType > > vectors_;
}; // class ContainerBasedDefault


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH
