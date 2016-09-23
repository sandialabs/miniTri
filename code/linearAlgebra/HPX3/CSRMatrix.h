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
// File:      CSRmatrix.h                                                   //
// Project:   miniTri                                                       //   
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Header file for CSR matrix class.                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef CSRMATRIX_H
#define CSRMATRIX_H

typedef enum {UNDEFINED,ADJACENCY,LOWERTRI,UPPERTRI,INCIDENCE} matrixtype;

#include <list>
#include <vector>
#include <map>

#include <hpx/hpx.hpp>
#include <boost/shared_array.hpp>

class Vector;

//////////////////////////////////////////////////////////////////////////////
// Compressed Sparse Row storage format Matrix
//////////////////////////////////////////////////////////////////////////////
class CSRMat 
{

 private:
  matrixtype type;
  int m;   //number of rows
  int n;   //number of cols
  //int nnz; //number of nonzeros

  boost::shared_array<int> nnzInRow;                             // nnz in each row
  boost::shared_array<boost::shared_array<int> > cols;           //columns of nonzeros
  boost::shared_array<boost::shared_array<int> > vals; //values of nonzeros
  boost::shared_array<boost::shared_array<int> > vals2; //values of nonzeros


  std::string mName;

  int mBlockSize;

 public:

  //Array of futures for matrix computation 
  //Maybe move to private eventually
  std::vector<hpx::shared_future<int> > mOps;




  //////////////////////////////////////////////////////////////////////////
  // default constructor -- builds empty matrix
  //////////////////////////////////////////////////////////////////////////
  CSRMat() 
    :type(UNDEFINED),m(0),n(0),nnzInRow(),cols(),vals(),vals2(),mName(""),mBlockSize(1)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Constructor that accepts matrix type as an argument
  //////////////////////////////////////////////////////////////////////////
  CSRMat(matrixtype _type,std::string str="",int bs=1) 
    :type(_type),m(0),n(0),nnzInRow(),cols(),vals(),vals2(),mName(str),mBlockSize(bs)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // constructor -- allocates memory for CSR sparse matrix
  //////////////////////////////////////////////////////////////////////////
   CSRMat(int _m, int _n, std::string str="",int bs=1, bool allocateVals2=false)
    :type(UNDEFINED),m(_m),n(_n),
    nnzInRow(new int[m]),cols(new boost::shared_array<int> [m]),
    vals(new boost::shared_array<int> [m]),
    mName(str),mBlockSize(bs)
  {

    if(allocateVals2==true)
    {
      vals2 = boost::shared_array<boost::shared_array<int> >(new boost::shared_array<int> [m]);
    }
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Copy constructor -- No data reallocation, just copying of smart pointers
  //////////////////////////////////////////////////////////////////////////
  CSRMat(const CSRMat &obj)
    //    :type(obj.type),m(obj.m),n(obj.n),nnz(obj.nnz),
    :type(obj.type),m(obj.m),n(obj.n),
    nnzInRow(obj.nnzInRow),cols(obj.cols), vals(obj.vals),vals2(obj.vals2),
    mName(obj.mName),mBlockSize(obj.mBlockSize)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // destructor -- deletes matrix
  //////////////////////////////////////////////////////////////////////////
  ~CSRMat()
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // additional functions and prototypes for additional functions
  // defined in CSRmatrix.cc
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // print -- prints matrix to stdio
  //////////////////////////////////////////////////////////////////
  void print() const;
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // Sums matrix elements
  //////////////////////////////////////////////////////////////////
  std::list<int> getSumElements() const;
  void getRowNNZs(std::vector<int> &rNNZs) const;
  void getColNNZs(std::vector<int> &cNNZs) const;
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // additional accessors and accessor prototypes
  //////////////////////////////////////////////////////////////////
  // returns the number of rows
  int getM() const { return m;};

  // returns the number of cols
  int getN() const { return n;};

  // returns the number of cols
  // Perhaps better way to do this, if nnz could be incremented by
  // each thread atomically
  int getNNZ() const 
  {
    int nnz=0;
    for(int i=0; i< m; i++)
    {
      nnz+=nnzInRow[i];
    } 
    return nnz;
  };

  // returns NNZ in row rnum
  inline int getNNZInRow(int rnum) const {return nnzInRow[rnum];};

  // returns column # for nonzero in row rowi at index nzindx
  inline int getCol(int rowi, int nzindx) const {return cols[rowi][nzindx];};


  // returns value for nonzero at inddex nzindx
  inline const int & getVal(int rowi, int nzindx) const 
    {return vals[rowi][nzindx];};

  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // y = this * 1
  //////////////////////////////////////////////////////////////////
  void SpMV1(bool trans, Vector &y);
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // level 3 basic linear algebra subroutines
  //////////////////////////////////////////////////////////////////
  void matmat(const CSRMat &A, const CSRMat &B);
  //////////////////////////////////////////////////////////////////

  void computeKCounts(const Vector &vTriDegrees,
		      const Vector &eTriDegrees, const std::map<int,std::map<int,int> > & edgeInds,
		      std::vector<int> &kCounts);

  void readMMMatrix(const char *fname);

  void createTriMatrix(const CSRMat &matrix, matrixtype mtype);
  void createIncidentMatrix(const CSRMat &matrix, std::map<int,std::map<int,int> > & eIndices);

  void permute();

  //////////////////////////////////////////////////////////////////////////

};
//////////////////////////////////////////////////////////////////////////////

#endif
