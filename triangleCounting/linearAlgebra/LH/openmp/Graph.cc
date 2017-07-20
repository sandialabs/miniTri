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
// Project:   miniTri                                                       //
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Source file for graph class.                                //
//////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sys/time.h>

#include "Graph.h"
#include "mmUtil.h"
#include "mmio.h"


//////////////////////////////////////////////////////////////////////////////
// Count triangles in graph
//////////////////////////////////////////////////////////////////////////////
void Graph::countTriangles()
{
  struct timeval t1, t2;
  double eTime;

  std::cout << "************************************************************"
            << "**********" << std::endl;
  std::cout << "Counting triangles ....." << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;

  ///////////////////////////////////////////////////////////////////////
  // Form B
  ///////////////////////////////////////////////////////////////////////
  std::cout << "--------------------" << std::endl;
  std::cout << "Creating incidence matrix B...";

  gettimeofday(&t1, NULL);

  CSRMat B(INCIDENCE,mBlockSize);
  B.createIncidenceMatrix(mMatrix,mEdgeIndices);

  gettimeofday(&t2, NULL);

  std::cout << " done" <<std::endl;

  //mMatrix.print();
  //B.print();

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
  std::cout << "TIME - Time to create B: " << eTime << std::endl;

  std::cout << "--------------------" << std::endl;
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // Create lower triangle matrix                                        
  ///////////////////////////////////////////////////////////////////////
  // Should eventually move this to constructor, build incidence from this
  CSRMat L(LOWERTRI);
  L.createTriMatrix(mMatrix, LOWERTRI);
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // C = L*B
  ///////////////////////////////////////////////////////////////////////
  std::cout << "--------------------" << std::endl;


  std::shared_ptr<CSRMat> C(new CSRMat(mMatrix.getM(),B.getN(),mBlockSize));

  std::cout << "C = L*B: " << std::endl;

  gettimeofday(&t1, NULL);
  C->matmat(L,B);
  gettimeofday(&t2, NULL);

  //C.print();

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
  std::cout << "TIME - Time to compute C = L*B: " << eTime << std::endl;

  std::cout << "--------------------" << std::endl;

  std::cout << "C NNZ: " << C->getNNZ() << std::endl;

  mNumTriangles = C->getNNZ();
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // Save triangle information                                           
  ///////////////////////////////////////////////////////////////////////
  mTriMat = C;
  ///////////////////////////////////////////////////////////////////////

  std::cout << "************************************************************"
            << "**********" << std::endl;
  std::cout << "Finished triangle counting" << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;
}
//////////////////////////////////////////////////////////////////////////////


