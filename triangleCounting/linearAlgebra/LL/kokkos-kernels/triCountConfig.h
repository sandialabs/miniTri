#ifndef _TRICOUNT_CONFIG_H 
#define _TRICOUNT_CONFIG_H 

#include  "KokkosSparse_CrsMatrix.hpp"

typedef Kokkos::OpenMP myExecSpace;

typedef int ordinal_t;
//typedef long long ordinal_t;
typedef double scalar_t;
typedef int size_type;
//typedef long long size_type;

typedef typename KokkosSparse::CrsMatrix<scalar_t, ordinal_t, myExecSpace, void, size_type > crsMat_t;
typedef typename crsMat_t::StaticCrsGraphType graph_t;
typedef typename crsMat_t::row_map_type::non_const_type row_map_view_t;
typedef typename crsMat_t::index_type::non_const_type   cols_view_t;
typedef typename crsMat_t::values_type::non_const_type values_view_t;

#endif
