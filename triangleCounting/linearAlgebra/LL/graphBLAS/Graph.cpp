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

//#include "GraphBLAS.h"

extern "C"{
#include "GraphBLAS.h"
}


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
  //  struct timeval t1, t2, t3, t4;
  double eTime1;

  //size_t numTriangles;
  int64_t numTriangles=0;


  std::cout << "************************************************************"
            << "**********" << std::endl;
  std::cout << "Counting triangles ....." << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;


  // Initialize Graph BLAS, should be moved to main?
  GrB_init(GrB_NONBLOCKING);


  //////////////////////////////////////////////////////////////////
  // typedefs
  //////////////////////////////////////////////////////////////////
  //  typedef typename crsMat_t::StaticCrsGraphType::row_map_type::non_const_type lno_view_t;
  typedef typename crsMat_t::StaticCrsGraphType::entries_type::non_const_type lno_nnz_view_t;
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Determine order of vertices/rows when forming L
  //////////////////////////////////////////////////////////////////
  const ordinal_t m = mAdjMatrixA.numRows();
  Kokkos::Impl::Timer timer1;

  std::vector<ordinal_t> newIndices;

  newIndices.resize(m);


  // Determine order of vertices/rows when forming L
  KokkosKernels::Impl::
    kk_sort_by_row_size <size_type, ordinal_t,myExecSpace>(m,mAdjMatrixA.graph.row_map.data(),
  							   newIndices.data());
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Create "lower triangle" matrix in graphBLAS format
  //////////////////////////////////////////////////////////////////
  size_type nnz = mAdjMatrixA.nnz();
  size_type nnzL = nnz/2;

  /////////////////////////////////////////////////////////////
  // Create tuple arrays i,j,v
  /////////////////////////////////////////////////////////////
  GrB_Index *i = (GrB_Index*) malloc (nnzL * sizeof (int64_t));
  GrB_Index *j = (GrB_Index*) malloc (nnzL * sizeof (int64_t));

  // Should this be int or something else?
  uint32_t *v = (uint32_t *) malloc (nnzL * sizeof (uint32_t));
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Copy data from KK structure to Graph BLAS arrays
  /////////////////////////////////////////////////////////////
  graph_t & graphA = mAdjMatrixA.graph; 

   auto rowmapA = graphA.row_map;
  //lno_view_t rowmapL = graphL.row_map;
  lno_nnz_view_t colIndsA = graphA.entries;

  size_type nzIndx=0;
  for(ordinal_t rownum=0; rownum<m; rownum++)
  {
    for(size_type indx=rowmapA[rownum]; indx<rowmapA[rownum+1]; indx++)
    {
      if(newIndices[colIndsA[indx]] < newIndices[rownum])
      {
        i[nzIndx] = rownum;
        j[nzIndx] = colIndsA[indx];
        v[nzIndx] = 1;
        nzIndx++;
      }
    }
  }
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Build Graph BLAS matrix from tuples
  /////////////////////////////////////////////////////////////
  GrB_Matrix gbL = NULL;

  GrB_Matrix_new (&gbL, GrB_UINT32, m, m);
  //For some reason this macro didn't seem to work
  //GrB_Matrix_build (C, i, j, v, nnz, GrB_PLUS_UINT32);
  GrB_Matrix_build_UINT32 (gbL, i, j, v, nnzL, GrB_PLUS_UINT32);
  /////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // C = L*L mask L using Graph BLAS
  ///////////////////////////////////////////////////////////////////////
  Kokkos::Impl::Timer timer2;
  GrB_Matrix gbC = NULL;
  //  OK (GrB_Matrix_nrows (&n, L)) ;
  GrB_Matrix_new (&gbC, GrB_UINT32, m, m);

  // gbC=gbL*gbL mask gbL                                                                                    
  GrB_mxm (gbC, gbL, NULL, GrB_PLUS_TIMES_UINT32, gbL, gbL, NULL);

  ///////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////
  // Sum up triangles
  ///////////////////////////////////////////////////////////////////////
  //GrB_reduce (&numTriangles, NULL, GrB_PLUS_INT64_MONOID, gbC, NULL);
  GrB_Matrix_reduce_INT64(&numTriangles, NULL, GrB_PLUS_INT64_MONOID, gbC, NULL);
  ///////////////////////////////////////////////////////////////////////

  myExecSpace::fence();
  eTime1 = timer1.seconds();


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

  std::cout << "Time to count triangles: " << eTime1 << std::endl;
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // free arrays
  //////////////////////////////////////////////////////////////////
  free(i);
  free(j);
  free(v);
  //////////////////////////////////////////////////////////////////



  GrB_finalize();


  return eTime1;
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


