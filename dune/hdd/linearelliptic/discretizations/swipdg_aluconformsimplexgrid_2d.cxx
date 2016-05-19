//#ifdef ENABLE_MPI
//# undef ENABLE_MPI
//#endif
//#define ENABLE_MPI 0

//#ifdef ENABLE_PARMETIS
//# undef ENABLE_PARMETIS
//#endif
//#define ENABLE_PARMETIS 0

//#include "config.h"

//#if HAVE_ALUGRID

//# include <dune/grid/alugrid.hh>

//# include "swipdg.hxx"

//namespace Dune {
//namespace HDD {
//namespace LinearElliptic {
//namespace Discretizations {


//template class SWIPDG< ALUGrid< 2, 2, simplex, conforming >,
//                       Stuff::Grid::ChooseLayer::leaf,
//                       double,
//                       1,
//                       1,
//                       GDT::ChooseSpaceBackend::fem,
//                       Stuff::LA::ChooseBackend::istl_sparse >;

//# if HAVE_DUNE_GRID_MULTISCALE

//template class SWIPDG< ALUGrid< 2, 2, simplex, conforming >,
//                       Stuff::Grid::ChooseLayer::local,
//                       double,
//                       1,
//                       1,
//                       GDT::ChooseSpaceBackend::fem,
//                       Stuff::LA::ChooseBackend::istl_sparse >;

//# endif // HAVE_DUNE_GRID_MULTISCALE


//} // namespace Discretizations
//} // namespace LinearElliptic
//} // namespace HDD
//} // namespace Dune

//#endif // HAVE_ALUGRID
