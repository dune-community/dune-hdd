// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH
#define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH

#ifndef DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BASE_DISABLE_CACHING
# define DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BASE_DISABLE_CACHING 0
#endif

#include <map>
#include <algorithm>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/pymor/operators/base.hh>
#include <dune/pymor/operators/affine.hh>
#include <dune/pymor/functionals/affine.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/products/boundaryl2.hh>
#include <dune/gdt/products/elliptic.hh>
#include <dune/gdt/products/h1.hh>
#include <dune/gdt/products/l2.hh>

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

  CachedDefault(TestSpaceType test_spc,
                AnsatzSpaceType ansatz_spc,
                Stuff::Common::Configuration bnd_inf_cfg,
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

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return test_space_;
  }

  const GridViewType& grid_view() const
  {
    return test_space_.grid_view();
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
    return VectorType(ansatz_space_.mapper().size());
  }

  void visualize(const VectorType& vector,
                 const std::string filename,
                 const std::string name,
                 const bool add_dirichlet = true,
                 Pymor::Parameter mu = Pymor::Parameter()) const
  {
    const auto vectors = this->available_vectors();
    if (add_dirichlet && std::find(vectors.begin(), vectors.end(), "dirichlet") != vectors.end()) {
      VectorType tmp = vector.copy();
      const auto dirichlet_vector = this->get_vector("dirichlet");
      if (dirichlet_vector.parametric()) {
        const Pymor::Parameter mu_dirichlet = this->map_parameter(mu, "dirichlet");
        if (mu_dirichlet.type() != dirichlet_vector.parameter_type())
          DUNE_THROW(Pymor::Exceptions::wrong_parameter_type,
                     mu_dirichlet.type() << " vs. " << dirichlet_vector.parameter_type());
        tmp += dirichlet_vector.freeze_parameter(mu);
      } else {
        tmp += *(dirichlet_vector.affine_part());
      }
      GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType >(ansatz_space_, tmp, name).visualize(filename);
    } else {
      GDT::ConstDiscreteFunction< AnsatzSpaceType, VectorType >(ansatz_space_, vector, name).visualize(filename);
    }
  } // ... visualize(...)

  void visualize(const VectorType& vector,
                 const std::string filename,
                 const std::string name,
                 Pymor::Parameter mu) const
  {
    visualize(vector, filename, name, true, mu);
  }

  using BaseType::solve;

  void solve(const DSC::Configuration options, VectorType& vector, const Pymor::Parameter mu = Pymor::Parameter()) const
  {
    auto logger = DSC::TimedLogger().get(static_id());
#if !DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BASE_DISABLE_CACHING
    const auto options_in_cache = cache_.find(options);
    bool exists = (options_in_cache != cache_.end());
    typename std::map< Pymor::Parameter, std::shared_ptr< VectorType > >::const_iterator options_and_mu_in_cache;
    if (exists) {
      options_and_mu_in_cache = options_in_cache->second.find(mu);
      exists = (options_and_mu_in_cache != options_in_cache->second.end());
    }
    if (!exists) {
#endif // !DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BASE_DISABLE_CACHING
      logger.info() << "solving";
      if (options.has_key("type"))
        logger.info() << " with '" << options.get< std::string >("type") << "'";
      if (!mu.empty())
        logger.info() << " for mu = " << mu;
      logger.info() << "... " << std::endl;
      uncached_solve(options, vector, mu);
#if !DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BASE_DISABLE_CACHING
      cache_[options][mu] = Stuff::Common::make_unique< VectorType >(vector.copy());
    } else {
      logger.info() << "retrieving solution ";
      if (!mu.empty())
        logger.info() << "for mu = " << mu << " ";
      logger.info() << "from cache... " << std::endl;
      const auto& result = *(options_and_mu_in_cache->second);
      vector = result;
    }
#endif // !DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_BASE_DISABLE_CACHING
  } // ... solve(...)

  void uncached_solve(const DSC::Configuration options, VectorType& vector, const Pymor::Parameter mu) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().uncached_solve(options, vector, mu));
  }

protected:
  const TestSpaceType test_space_;
  const AnsatzSpaceType ansatz_space_;
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
  using typename BaseType::PatternType;
  typedef typename VectorType::ScalarType RangeFieldType;

protected:
  typedef Pymor::LA::AffinelyDecomposedContainer< MatrixType > AffinelyDecomposedMatrixType;
  typedef Pymor::LA::AffinelyDecomposedContainer< VectorType > AffinelyDecomposedVectorType;
  typedef Stuff::LA::Solver< MatrixType > SolverType;

public:
  static std::string static_id() { return "hdd.linearelliptic.discretizations.containerbased"; }

  ContainerBasedDefault(TestSpaceType test_spc,
                        AnsatzSpaceType ansatz_spc,
                        const Stuff::Common::Configuration& bnd_inf_cfg,
                        const ProblemType& prb)
    : BaseType(test_spc, ansatz_spc, bnd_inf_cfg, prb)
    , container_based_initialized_(false)
    , purely_neumann_(false)
    , matrix_(std::make_shared< AffinelyDecomposedMatrixType >())
    , rhs_(std::make_shared< AffinelyDecomposedVectorType >())
    , pattern_(nullptr)
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

  const PatternType& pattern() const
  {
    return *pattern_;
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

  std::vector< std::string > solver_types() const
  {
    return SolverType::types();
  }

  DSC::Configuration solver_options(const std::string type = "") const
  {
    return SolverType::options(type);
  }

  /**
   * \brief solves for u_0
   */
  void uncached_solve(const DSC::Configuration options,
                      VectorType& vector,
                      const Pymor::Parameter mu = Pymor::Parameter()) const
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
      SolverType(tmp_system_matrix).apply(tmp_rhs, vector, options);
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
        frozenOperator.apply_inverse(*rhs_vector, vector, options);
      } else {
        const auto nonparametricOperator = lhsOperator.affine_part();
        nonparametricOperator.apply_inverse(*rhs_vector, vector, options);
      }
    }
  } // ... uncached_solve(...)

protected:
  template< class SetConstraints, class ClearConstraints >
  std::shared_ptr< AffinelyDecomposedMatrixType >
  make_zero_dirichlet_product(const std::shared_ptr< AffinelyDecomposedMatrixType >& prod,
                              const SetConstraints& set,
                              const ClearConstraints& clear)
  {
    auto ret = std::make_shared<AffinelyDecomposedMatrixType>(prod->copy());
    for (ssize_t qq = 0; qq < ret->num_components(); ++qq) {
      // clear rows
      clear.apply(*ret->component(qq));
      // clear columns
      for (const auto& column : clear.dirichlet_DoFs())
        ret->component(qq)->unit_col(column);
    }
    if (!ret->has_affine_part()) {
      // we can assume that there is at least one component
      assert(ret->num_components() > 0);
      ret->register_affine_part(new MatrixType(ret->component(0)->copy()));
      *ret->affine_part() *= 0;
    }
    // set rows
    set.apply(*ret->affine_part());
    // set columns
    for (const auto& column : set.dirichlet_DoFs())
      ret->affine_part()->unit_col(column);
    return ret;
  } // ... make_zero_dirichlet_product(...)

  void assemble_products(const std::vector< std::string > only_these_products, const size_t over_integrate = 2)
  {
    using namespace GDT;

    if (only_these_products.size() == 0)
      return;

    GDT::SystemAssembler< TestSpaceType, GridViewType, AnsatzSpaceType > system_assembler(this->test_space(),
                                                                                          this->ansatz_space(),
                                                                                          this->grid_view());
    // L2
    typedef Products::L2Assemblable< MatrixType, TestSpaceType, GridViewType, AnsatzSpaceType > L2ProductType;
    std::unique_ptr< L2ProductType > l2_product;
    auto l2_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
    if (   std::find(only_these_products.begin(), only_these_products.end(), "l2") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "l2_0") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "h1") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "h1_0") != only_these_products.end()) {
      l2_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                              this->ansatz_space().mapper().size(),
                                              *pattern_);
      l2_product = DSC::make_unique< L2ProductType >(*l2_product_matrix->affine_part(),
                                                     this->test_space(),
                                                     this->grid_view(),
                                                     this->ansatz_space(),
                                                     over_integrate);
      system_assembler.add(*l2_product);
    }
    // H1 semi
    typedef Products::H1SemiAssemblable< MatrixType, TestSpaceType, GridViewType, AnsatzSpaceType >
        SemiH1ProductType;
    std::unique_ptr< SemiH1ProductType > semi_h1_product;
    auto semi_h1_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
    if (   std::find(only_these_products.begin(), only_these_products.end(), "h1_semi") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "h1_semi_0") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "h1") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "h1_0") != only_these_products.end()) {
      semi_h1_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                   this->ansatz_space().mapper().size(),
                                                   *pattern_);
      semi_h1_product = DSC::make_unique< SemiH1ProductType >(*(semi_h1_product_matrix->affine_part()),
                                                              this->test_space(),
                                                              this->grid_view(),
                                                              this->ansatz_space(),
                                                              over_integrate);
      system_assembler.add(*semi_h1_product);
    }
    // elliptic
    const auto& diffusion_factor = *this->problem().diffusion_factor();
    const auto& diffusion_tensor = *this->problem().diffusion_tensor();
    assert(!diffusion_tensor.parametric());
    typedef Products::EllipticAssemblable< MatrixType, typename ProblemType::DiffusionFactorType::NonparametricType,
                                           TestSpaceType, GridViewType, AnsatzSpaceType, RangeFieldType,
                                           typename ProblemType::DiffusionTensorType::NonparametricType > EllipticProductType;
    std::vector< std::unique_ptr< EllipticProductType > > elliptic_products;
    auto elliptic_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
    if (   std::find(only_these_products.begin(), only_these_products.end(), "elliptic") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "elliptic_0") != only_these_products.end()) {
      for (DUNE_STUFF_SSIZE_T qq = 0; qq < diffusion_factor.num_components(); ++qq) {
        const auto id = elliptic_product_matrix->register_component(diffusion_factor.coefficient(qq),
                                                                    this->test_space().mapper().size(),
                                                                    this->ansatz_space().mapper().size(),
                                                                    *pattern_);
        elliptic_products.emplace_back(new EllipticProductType(*elliptic_product_matrix->component(id),
                                                               this->test_space(),
                                                               this->grid_view(),
                                                               this->ansatz_space(),
                                                               *diffusion_factor.component(qq),
                                                               *diffusion_tensor.affine_part(),
                                                               over_integrate));
      }
      if (diffusion_factor.has_affine_part()) {
        elliptic_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                      this->ansatz_space().mapper().size(),
                                                      *pattern_);
        elliptic_products.emplace_back(new EllipticProductType(*elliptic_product_matrix->affine_part(),
                                                               this->test_space(),
                                                               this->grid_view(),
                                                               this->ansatz_space(),
                                                               *diffusion_factor.affine_part(),
                                                               *diffusion_tensor.affine_part(),
                                                               over_integrate));
      }
      for (auto& product : elliptic_products)
        system_assembler.add(*product);
    }
    // boundary L2
    typedef Products::BoundaryL2Assemblable< MatrixType, TestSpaceType, GridViewType, AnsatzSpaceType >
        BoundaryL2ProductType;
    std::unique_ptr< BoundaryL2ProductType > boundary_l2_product;
    auto boundary_l2_product_matrix = std::make_shared< AffinelyDecomposedMatrixType >();
    if (   std::find(only_these_products.begin(), only_these_products.end(), "boundary_l2") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "boundary_l2_0") != only_these_products.end()) {
      boundary_l2_product_matrix->register_affine_part(this->test_space().mapper().size(),
                                                       this->ansatz_space().mapper().size(),
                                                       *pattern_);
      boundary_l2_product = DSC::make_unique< BoundaryL2ProductType >(*(boundary_l2_product_matrix->affine_part()),
                                                                      this->test_space(),
                                                                      this->grid_view(),
                                                                      this->ansatz_space(),
                                                                      over_integrate);
      system_assembler.add(*boundary_l2_product);
    }
    // walk the grid
    system_assembler.assemble();
    // register the products
    for (auto&& key_value_pair : {std::make_pair("l2", l2_product_matrix),
                                  std::make_pair("h1_semi", semi_h1_product_matrix),
                                  std::make_pair("elliptic", elliptic_product_matrix)})  {
      auto key = std::string(key_value_pair.first);
      auto value = key_value_pair.second;
      if (   std::find(only_these_products.begin(), only_these_products.end(), key) != only_these_products.end()
          || std::find(only_these_products.begin(), only_these_products.end(), key + "_0") != only_these_products.end())
        products_.insert(std::make_pair(key, value));
    }
    // energy product is the system matrix
    if (   std::find(only_these_products.begin(), only_these_products.end(), "energy") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "energy_") != only_these_products.end())
      products_.insert(std::make_pair("energy", std::make_shared<AffinelyDecomposedMatrixType>(matrix_->copy())));
    // h1 product is l2 product + h1_semi product
    if (   std::find(only_these_products.begin(), only_these_products.end(), "h1") != only_these_products.end()
        || std::find(only_these_products.begin(), only_these_products.end(), "h1_0") != only_these_products.end()) {
      auto h1_product_matrix
          = std::make_shared< AffinelyDecomposedMatrixType >(new MatrixType(l2_product_matrix->affine_part()->backend()));
      h1_product_matrix->affine_part()->axpy(1.0, *semi_h1_product_matrix->affine_part());
      products_.insert(std::make_pair("h1", h1_product_matrix));
    }
  } // ... assemble_products(...)

  template< class I >
  void assemble_products(const std::vector< std::string > only_these_products,
                         const GDT::Spaces::DirichletConstraints< I >& clear_and_set_dirichlet_rows,
                         const GDT::Spaces::DirichletConstraints< I >& clear_dirichlet_rows,
                         const size_t over_integrate = 2)
  {
    using namespace GDT;

    if (only_these_products.size() == 0)
      return;

    assemble_products(only_these_products, over_integrate);

    for (const std::string& prod_0 : only_these_products) {
      if (prod_0.size() > 2 && prod_0.substr(prod_0.size() - 2, 2) == "_0") {
        const auto prod = prod_0.substr(0, prod_0.size() - 2);
        products_.insert(std::make_pair(prod_0, make_zero_dirichlet_product(products_.at(prod),
                                                                            clear_and_set_dirichlet_rows,
                                                                            clear_dirichlet_rows)));
      }
    }
  } // assemble_product(...)

  void finalize_init(const bool prune)
  {
    if (!container_based_initialized_) {
      if (prune) {
        matrix_ = std::make_shared< AffinelyDecomposedMatrixType >(matrix_->pruned());
        for (auto& element : products_)
          element.second = std::make_shared< AffinelyDecomposedMatrixType >(element.second->pruned());
      }
      container_based_initialized_ = true;
    }
  } // ... finalize_init(...)

  void assert_everything_is_ready() const
  {
    if (!container_based_initialized_)
      DUNE_THROW(Stuff::Exceptions::you_are_using_this_wrong,
                 "The implemented discretization has to fill 'matrix_', 'rhs_' and 'pattern_' during init() and call "
                 << "finalize_init()!\n"
                 << "The user has to call init() before calling any other method!");
  } // ... assert_everything_is_ready()

  bool container_based_initialized_;
  bool purely_neumann_;
  std::shared_ptr< AffinelyDecomposedMatrixType > matrix_;
  std::shared_ptr< AffinelyDecomposedVectorType > rhs_;
  std::shared_ptr< PatternType > pattern_;
  mutable std::map< std::string, std::shared_ptr< AffinelyDecomposedMatrixType > > products_;
  mutable std::map< std::string, std::shared_ptr< AffinelyDecomposedVectorType > > vectors_;
}; // class ContainerBasedDefault


} // namespace Discretizations
} // namespace LinearElliptic
} // namespace HDD
} // namespace Dune

#endif // DUNE_HDD_LINEARELLIPTIC_DISCRETIZATIONS_DEFAULT_HH
