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
// File:      Graph.h                                                       //
// Project:   miniTri: triangleCounting                                     //   
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Header file for graph class.                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef GRAPH_H
#define GRAPH_H

#include "triCountConfig.h"
#include "KokkosSparse_CrsMatrix.hpp"

void readMMMatrixA(const char* fname, crsMat_t &outMat);

//////////////////////////////////////////////////////////////////////////////
// Graph class
//////////////////////////////////////////////////////////////////////////////
class Graph 
{

 private:

  ordinal_t mNumVerts;
  ordinal_t mNumEdges;

  crsMat_t mAdjMatrixA;

 public:
  //////////////////////////////////////////////////////////////////////////
  // default constructor -- builds empty graph
  //////////////////////////////////////////////////////////////////////////
  Graph() 
    :mNumVerts(0),mAdjMatrixA()
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Constructor that accepts matrix type as an argument
  //////////////////////////////////////////////////////////////////////////
  Graph(std::string fname1) 
    :mAdjMatrixA()
  {
    struct timeval t1, t2;

    gettimeofday(&t1, NULL);

    ////////////////////////////////////////////////////////////////////
    // Reading input file
    ////////////////////////////////////////////////////////////////////

    std::cout << "Reading input file..." << std::endl;
    readMMMatrixA(fname1.c_str(),mAdjMatrixA);

    std::cout << "Finished Reading input file." << std::endl;
    ////////////////////////////////////////////////////////////////////
    gettimeofday(&t2, NULL);


    double eTime;
    eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
    std::cout << "Time to read input file: " << eTime << std::endl;

  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // destructor -- deletes matrix
  //////////////////////////////////////////////////////////////////////////
  ~Graph()
  {
  };
  //////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////
  // Count triangles
  //////////////////////////////////////////////////////////////////////////
  double triangleCount();
  //////////////////////////////////////////////////////////////////////////

};
//////////////////////////////////////////////////////////////////////////////

#endif
