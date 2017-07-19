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
//#include "util.h"
#include "mmio.h"

//////////////////////////////////////////////////////////////////////////////
// Enumerate triangles in graph
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
  ///////////////////////////////////////////////////////////////////////
  std::cout << "--------------------" << std::endl;
  std::cout << "Creating L,U matrices ...";

  gettimeofday(&t1, NULL);

  CSRMat L(LOWERTRI);
  L.createTriMatrix(mMatrix, LOWERTRI);

  CSRMat U(UPPERTRI);
  U.createTriMatrix(mMatrix, UPPERTRI);

  gettimeofday(&t2, NULL);

  std::cout << " done" <<std::endl;

  //L.print();
  //U.print();

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
  std::cout << "TIME - Time to create L, U  matrices: " << eTime << std::endl;

  std::cout << "--------------------" << std::endl;
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // B = L*U
  ///////////////////////////////////////////////////////////////////////
  std::cout << "--------------------" << std::endl;

  CSRMat B(L.getM(),U.getN());

  std::cout << "B = L*U: " << std::endl;

  gettimeofday(&t1, NULL);
  B.matmat(L,U);
  gettimeofday(&t2, NULL);

  //  B.print();

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
  std::cout << "TIME - Time to compute B = L*U: " << eTime << std::endl;

  std::cout << "--------------------" << std::endl;

  std::cout << "NNZ: " << B.getNNZ() << std::endl;

  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // B = B .* L
  ///////////////////////////////////////////////////////////////////////
  std::cout << "--------------------" << std::endl;
  std::cout << "B = B .* L: " << std::endl;

  gettimeofday(&t1, NULL);
  B.EWMult(L);
  gettimeofday(&t2, NULL);

  //  B.print();

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
  std::cout << "TIME - Time to compute B = B .* L: " << eTime << std::endl;
  std::cout << "--------------------" << std::endl;
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // Counting triangles
  ///////////////////////////////////////////////////////////////////////
  std::cout << "--------------------" << std::endl;

  gettimeofday(&t1, NULL);
  mTriangles = B.getSumElements();
  gettimeofday(&t2, NULL);

  assert(mTriangles.size() % 3 == 0);
  mNumTriangles = mTriangles.size() / 3;

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
  std::cout << "TIME - Time to sum up triangles: " << eTime << std::endl;

  std::cout << "--------------------" << std::endl;
  ///////////////////////////////////////////////////////////////////////


  std::cout << "************************************************************"
            << "**********" << std::endl;
  std::cout << "Finished triangle counting" << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;
}
//////////////////////////////////////////////////////////////////////////////

