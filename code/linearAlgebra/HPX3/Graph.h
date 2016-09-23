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
// Project:   miniTri                                                       //
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Header file for Graph class.                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef GRAPH_H
#define GRAPH_H

#include <list>
#include <vector>
#include <map>
#include <cmath>

//#include "mmio.h"
#include "CSRMatrix.h"
#include "Vector.h"

#include <boost/shared_ptr.hpp>

//////////////////////////////////////////////////////////////////////////////
// Graph class
//////////////////////////////////////////////////////////////////////////////
class Graph 
{

 private:
  std::string mFilename;

  int mNumVerts;
  int mNumEdges;
  CSRMat mMatrix;

  boost::shared_ptr<CSRMat> mTriMat; //Matrix that contains triangle info

  std::map<int,std::map<int,int> > mEdgeIndices;

  // Triangle degrees
  Vector mVTriDegrees;
  Vector mETriDegrees;

  // Size of granularity used by HPX calls
  int mBlockSize;

  // K-count frequency table
  std::vector<int> mKCounts;

 public:
  //////////////////////////////////////////////////////////////////////////
  // default constructor -- builds empty graph
  //////////////////////////////////////////////////////////////////////////
  Graph() 
    :mFilename("UNDEFINED"),mNumVerts(0),mMatrix(), mTriMat(),mBlockSize(1)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Constructor that accepts matrix type as an argument
  //////////////////////////////////////////////////////////////////////////
  Graph(std::string _fname, int bs=1) 
    :mFilename(_fname),mMatrix(ADJACENCY,"A",bs), mTriMat(),mBlockSize(bs)
  {
     mMatrix.readMMMatrix(mFilename.c_str());
     mNumVerts = mMatrix.getM();
     mNumEdges = mMatrix.getNNZ() /2;

     std::cout << "Num Vertices: " << mNumVerts << std::endl;
     std::cout << "Num Edges: " << mNumEdges << std::endl;

     int countSize = (int) sqrt(mNumVerts);
     if(countSize < 10)  // Ask Jon how to set this
     {
       countSize = 10;
     }
     mKCounts.resize(countSize);

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
  // Enumerate triangles
  //////////////////////////////////////////////////////////////////////////
  void triangleEnumerate();
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Calculate triangle degrees
  //////////////////////////////////////////////////////////////////////////
  void calculateTriangleDegrees();
  //////////////////////////////////////////////////////////////////////////

  void calculateKCounts();
  void printKCounts();

  void printTriangles() const;

  // This call introduces global barrier
  void printNumTriangles() const
  {
    hpx::wait_all(mTriMat->mOps);
    // Account for 3 instances of each triangle
    int numTris = mTriMat->getNNZ()/3;

    std::cout << "WARNING -- GLOBAL BARRIER --" << std::endl;
    std::cout << "   Number of triangles: " << numTris << std::endl;
  };



};
//////////////////////////////////////////////////////////////////////////////

#endif
