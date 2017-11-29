//@HEADER
// ************************************************************************
// 
//                        miniTri v. 1.0
//              Copyright (2016) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  Jon Berry (jberry@sandia.gov)
//                     Michael Wolf (mmwolf@sandia.gov)
// 
// ************************************************************************
//@HEADER

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File:      Graph.cc                                                      //
// Project:   miniTri: triangle counting                                    //
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Source file for graph class.                                //
//////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <sys/time.h>

#include "Graph.h"
#include "mmUtil.h"
#include "mmio.h"

#include <impl/Kokkos_Timer.hpp>

#include "KokkosGraph_Triangle.hpp"

//////////////////////////////////////////////////////////////////////////////
// Count triangles in graph
//////////////////////////////////////////////////////////////////////////////
double Graph::triangleCount()
{
  struct timeval t1, t2, t3, t4;
  double eTime;

  size_t numTriangles;


  std::cout << "************************************************************"
            << "**********" << std::endl;
  std::cout << "Counting triangles ....." << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;



  //////////////////////////////////////////////////////////////////
  // typedefs
  //////////////////////////////////////////////////////////////////
  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;
  typedef myExecSpace TempMemSpace;
  typedef myExecSpace PersistentMemSpace;

  typedef KokkosKernels::Experimental::KokkosKernelsHandle
    <lno_view_t,lno_nnz_view_t, lno_nnz_view_t, myExecSpace, TempMemSpace,PersistentMemSpace > KernelHandle;
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Set Parameters to tune kkmem algorithm
  //////////////////////////////////////////////////////////////////

  // Handle for L*L, mask L triangle enumeration algorithm
  KernelHandle kh;
  kh.create_spgemm_handle(KokkosSparse::SPGEMM_KK_TRIANGLE_LL);

  kh.set_dynamic_scheduling(true);


  // Sets thread chunk size
  int chunk_size=16;
  kh.set_team_work_size(chunk_size);

  // 2 steps of compression (true would be 1 step of compression)
  kh.get_spgemm_handle()->set_compression_steps(false);

  kh.set_verbose(true);


  //////////////////////////////////////////////////////
  // Other Kokkos Kernels options that were not used
  //////////////////////////////////////////////////////

  //int shmemsize =100000;
  //kh.set_shmem_size(shmemsize);
  //kh.set_suggested_team_size(team_size);

  // Tells SGEMM not to sort the lower triangular since we are going to sort this ourselves
  // kh.get_spgemm_handle()->set_sort_lower_triangular(false);

  //This is to calculate the number of operations
  //kh.get_spgemm_handle()->set_read_write_cost_calc(true);

  //////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Allocate memory for output structure
  //////////////////////////////////////////////////////////////////
  const ordinal_t m = mAdjMatrixA.numRows();

  Kokkos::View <size_t *,myExecSpace> matC_row;

  matC_row = Kokkos::View <size_t *,myExecSpace> ("Counts of triangles for each row", m);
  //////////////////////////////////////////////////////////////////

  //  crsMat_t adjMatrixL;

  // std::vector<ordinal_t> newIndices;

  // newIndices.resize(m);



  ///////////////////////////////////////////////////////////////////////
  // C = L*L mask L
  ///////////////////////////////////////////////////////////////////////
  Kokkos::Impl::Timer timer1;

  // For now comment these out, eventually want to explicitly do this
  //
  // //////////////////////////////////////////////////////////////////
  // // Creating lower triangular matrix.
  // //////////////////////////////////////////////////////////////////

  // // Determine order of vertices/rows when forming L
  // KokkosKernels::Experimental::Util::
  //   kk_sort_by_row_size <size_type, ordinal_t,myExecSpace>(m,mAdjMatrixA.graph.row_map.data(),
  // 							   newIndices.data());

  // // Create "lower triangular" portion of matrix
  // adjMatrixL = KokkosKernels::Experimental::Util::kk_get_lower_crs_matrix(mAdjMatrixA,newIndices.data());
  // //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Triangle counting
  // Counts triangles using linear algebra formulation for C = L*L with mask
  //
  // kh -- kernel handle
  // m -- Number of rows in L and C
  // graphL.row_map -- KokkosView of row map, matrix L
  // graphL.entries -- KokkosView of column indices, matrix L
  // Q kokkos Lambda -- run on each entry of resulting matrix
  //////////////////////////////////////////////////////////////////
  //  graph_t & graphL = adjMatrixL.graph;
  graph_t & graphA = mAdjMatrixA.graph;

  KokkosGraph::Experimental::triangle_generic(&kh, m, graphA.row_map, graphA.entries,
      KOKKOS_LAMBDA(const ordinal_t& row, const ordinal_t &col_set_index, const ordinal_t &col_set,  const ordinal_t &thread_id) 
      {
	matC_row(row) += KokkosKernels::Impl::pop_count(col_set);
      }
    );
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Sum up triangles
  //////////////////////////////////////////////////////////////////
  numTriangles = 0;
  KokkosKernels::Impl::
  kk_reduce_view<Kokkos::View <size_t *,myExecSpace>, myExecSpace>(m, matC_row, numTriangles);
  //////////////////////////////////////////////////////////////////

  myExecSpace::fence();
  eTime = timer1.seconds();


  std::cout << "************************************************************"
            << "**********" << std::endl;
  std::cout << "Finished triangle counting" << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;

  //////////////////////////////////////////////////////////////////
  // Output triangle counting info
  //////////////////////////////////////////////////////////////////
  ordinal_t numVerts = mAdjMatrixA.numRows();
  size_type numEdges = mAdjMatrixA.nnz()/2;

  std::cout << "|V| = " << numVerts << std::endl;
  std::cout << "|E| = " << numEdges << std::endl;
  std::cout << "|T| = " << numTriangles << std::endl;

  std::cout << "Time to count triangles: " << eTime << std::endl;
  //////////////////////////////////////////////////////////////////


  return eTime;
}
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Reads matrix market file
////////////////////////////////////////////////////////////////////////////////
void readMMMatrixA(const char *fname, crsMat_t &matrix)
{
   matrix = KokkosKernels::Impl::read_kokkos_crst_matrix<crsMat_t>(fname);
   matrix.values = values_view_t();
}
//////////////////////////////////////////////////////////////////////////////


