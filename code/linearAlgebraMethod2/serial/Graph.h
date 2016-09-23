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
//              Header file for graph class.                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef GRAPH_H
#define GRAPH_H

#include <list>
#include <vector>
#include <map>
#include <cmath>

#include "CSRmatrix.h"

//////////////////////////////////////////////////////////////////////////////
// Graph class
//////////////////////////////////////////////////////////////////////////////
class Graph 
{

 private:
  std::string mFilename;

  std::vector<int> mVertProp;

  int mNumVerts;
  CSRMat mMatrix;

  int mNumTriangles;
  std::list<int> mTriangles; //contains triangles in graph

  // Degree info                                                                                                                                                            
  std::map<int,int> mVDegrees;
  std::map<int,std::map<int,int> > mEDegrees;

  // K-count frequency table                                                                                                                                                
  std::vector<int> mKCounts;


 public:
  //////////////////////////////////////////////////////////////////////////
  // default constructor -- builds empty graph
  //////////////////////////////////////////////////////////////////////////
  Graph() 
    :mFilename("UNDEFINED"),mVertProp(0),mNumVerts(0),mMatrix(), mNumTriangles(0), mTriangles(0),
    mVDegrees(), mEDegrees(), mKCounts(0)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Constructor that accepts matrix type as an argument
  //////////////////////////////////////////////////////////////////////////
  Graph(std::string _fname) 
    :mFilename(_fname), mVertProp(0),mMatrix(), mNumTriangles(0), mTriangles(0),
    mVDegrees(), mEDegrees(), mKCounts(0)
  {
     mMatrix.readMMMatrix(mFilename.c_str());
     mNumVerts = mMatrix.getM();

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

  void printTriangles() const;

  int getNumTriangles() const {return mNumTriangles;};

  void calculateTriangleDegrees();
  void orderTriangles();
  void calculateKCounts();
  void printKCounts();

  void findMinTriDegrees(int v1, int v2, int v3, unsigned int &tvMin, unsigned int &teMin);

};
//////////////////////////////////////////////////////////////////////////////

#endif
