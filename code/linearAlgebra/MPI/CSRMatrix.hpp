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

typedef enum {UNDEFINED,LOWERTRI,UPPERTRI} matrixtype;

#include <list>
#include <vector>
#include <map>

#include <boost/shared_array.hpp>
#include <mpi.h>

#include "mmUtil.h"

class Vector;


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
struct CSRSubmat
{
  std::vector<std::vector<int> > cols;  //columns of nonzeros                                
  std::vector<std::vector<int> > vals;  //values of nonzeros                                 
  std::vector<std::vector<int> > vals2; //values of nonzeros                                 
};
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
// Compressed Sparse Row storage format Matrix
//////////////////////////////////////////////////////////////////////////////
class CSRMat 
{

 private:
  matrixtype type;
  int mGlobNumRows;   //number of global rows
  int mGlobNumCols;   //number of global cols
  int mGlobNNZ;       //number of global nonzeros
  int mLocNumRows;
  int mLocNNZ;
  int mStartRow;
  std::vector<int> mStartRowOnProc;

  std::vector<CSRSubmat> mSubmat;
  std::vector<int> mSubmatNumCols;
  std::vector<int> mSubmatStartCols;

  // MPI info
  MPI_Comm mComm;
  int mWorldSize;
  int mMyRank;


 public:
  //////////////////////////////////////////////////////////////////////////
  // default constructor -- builds empty matrix
  //////////////////////////////////////////////////////////////////////////
  CSRMat(MPI_Comm _comm=MPI_COMM_WORLD) 
    :type(UNDEFINED),mGlobNumRows(0),mGlobNumCols(0),mGlobNNZ(0),mLocNumRows(0),
     mLocNNZ(0), mComm(_comm)
  {
    MPI_Comm_size(mComm,&mWorldSize);
    MPI_Comm_rank(mComm,&mMyRank);

    mStartRowOnProc.resize(mWorldSize);
    mSubmat.resize(mWorldSize);
    mSubmatNumCols.resize(mWorldSize);
    mSubmatStartCols.resize(mWorldSize);
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Constructor that accepts matrix type as an argument
  //////////////////////////////////////////////////////////////////////////
  CSRMat(matrixtype _type,MPI_Comm _comm=MPI_COMM_WORLD) 
    :type(_type),mGlobNumRows(0),mGlobNumCols(0),mGlobNNZ(0),mLocNumRows(0), 
     mLocNNZ(0), mComm(_comm)
  {
    MPI_Comm_size(mComm,&mWorldSize);
    MPI_Comm_rank(mComm,&mMyRank);

    mStartRowOnProc.resize(mWorldSize);
    mSubmat.resize(mWorldSize);
    mSubmatNumCols.resize(mWorldSize);
    mSubmatStartCols.resize(mWorldSize);
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // constructor -- allocates memory for CSR sparse matrix
  //////////////////////////////////////////////////////////////////////////
  CSRMat(int _m, int _n, MPI_Comm _comm, bool allocateVals2=false)
    :type(UNDEFINED),mGlobNumRows(_m),mGlobNumCols(_n),mGlobNNZ(0),
     mComm(_comm)
  {
    MPI_Comm_size(mComm,&mWorldSize);
    MPI_Comm_rank(mComm,&mMyRank);

    partitionMatrix(mGlobNumRows,mWorldSize,mMyRank,mLocNumRows,mStartRow);

    mSubmat.resize(mWorldSize);
    for(int submatNum=0;submatNum<mWorldSize;submatNum++)
    {
      mSubmat[submatNum].cols.resize(mLocNumRows);
      mSubmat[submatNum].vals.resize(mLocNumRows);

      if(allocateVals2==true)
      {
        mSubmat[submatNum].vals.resize(mLocNumRows);
      }

    }

    mStartRowOnProc.resize(mWorldSize);
    mSubmatNumCols.resize(mWorldSize);
    mSubmatStartCols.resize(mWorldSize);

    //////////////////////////////////////////////////////////////////////
    // Sets values for mStartRowOnProc
    //////////////////////////////////////////////////////////////////////
    for(int rank=0; rank<mWorldSize; rank++)
    {
      int numrows;
      partitionMatrix(mGlobNumRows,mWorldSize,rank,
	  	      numrows,mStartRowOnProc[rank]);
    }
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Sets values for mSubmatStartCols,mSubmatNumCols                                        
    //////////////////////////////////////////////////////////////////////
    for(int rank=0; rank<mWorldSize; rank++)
    {
      partitionMatrix(mGlobNumCols,mWorldSize,rank,
	  	      mSubmatNumCols[rank],mSubmatStartCols[rank]);
    }
    //////////////////////////////////////////////////////////////////////

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
  //////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // additional accessors and accessor prototypes
  //////////////////////////////////////////////////////////////////

  // returns the global number of rows
  int getGlobNumRows() const { return mGlobNumRows;};

  // returns the global number of cols
  int getGlobNumCols() const { return mGlobNumCols;};

  // returns the global number of cols
  int getGlobNNZ() const { return mGlobNNZ;};

  // returns the local number of cols
  int getLocNumRows() const { return mLocNumRows;};

  // returns the start row
  int getStartRow() const { return mStartRow;};

  // returns local number of range elements
  int getLocNumRangeIDs() const 
  {
    return mSubmatNumCols[mMyRank];
  }

  // returns NNZ in row rnum
  inline int getNNZInRow(int submatNum, int rnum) const
  {return mSubmat[submatNum].cols[rnum].size();};

  // returns column # for nonzero in row rowi at index nzindx
  inline int getCol(int submatNum, int rowi, int nzindx) const
  {return mSubmat[submatNum].cols[rowi][nzindx];};

  inline const CSRSubmat & getSubMatrix(int submatNum) const
  { return mSubmat[submatNum]; }


  // returns true if row is on this process
  inline bool rowOnProc(int rowi) const
  {
    if(rowi >= mStartRow && rowi < mStartRow+mLocNumRows)
    {
      return true;
    }
    return false;
  }

  int  whichSubMatrix(int colID);
  int  rowOnWhichProc(int rowID);


  //////////////////////////////////////////////////////////////////

  void SpMV1(bool trans, Vector &y);

  //////////////////////////////////////////////////////////////////
  // level 3 basic linear algebra subroutines
  //////////////////////////////////////////////////////////////////
  void matmat(const CSRMat &A, const CSRMat &B);
  //////////////////////////////////////////////////////////////////

  void computeKCounts(const Vector &vTriDegrees, const Vector &eTriDegrees,
                      std::map<int,std::map<int,int> > & edgeInds,
                      std::vector<int> &kCounts);

  void readMMMatrix(const char* fname);
  void readBinMatrix(const char* fname);

  void createIncidenceMatrix(const CSRMat &matrix, std::map<int,std::map<int,int> > & eIndices);
  //////////////////////////////////////////////////////////////////////////

};
//////////////////////////////////////////////////////////////////////////////

#endif
