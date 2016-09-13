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
// File:      Vector.h                                                      //
// Project:   miniTri                                                       //   
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Header file for vector class.                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <mpi.h>

//////////////////////////////////////////////////////////////////////////////
// Compressed Sparse Row storage format Matrix
//////////////////////////////////////////////////////////////////////////////
class Vector
{

 private:

  std::vector<int> mElements;                             // elements in vector
  MPI_Comm mComm;
  int mWorldSize;
  int mMyRank;

 public:

  //////////////////////////////////////////////////////////////////////////
  // default constructor -- builds empty matrix
  //////////////////////////////////////////////////////////////////////////
  Vector(MPI_Comm _comm) 
    :mElements(),mComm(_comm)
  {
    MPI_Comm_size(mComm,&mWorldSize);
    MPI_Comm_rank(mComm,&mMyRank);
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // constructor -- allocates memory for CSR sparse matrix
  //////////////////////////////////////////////////////////////////////////
  Vector(int _m, MPI_Comm _comm)
    :mElements(_m), mComm(_comm)
  {
    MPI_Comm_size(mComm,&mWorldSize);
    MPI_Comm_rank(mComm,&mMyRank);
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Copy constructor -- No data reallocation, just copying of smart pointers
  //////////////////////////////////////////////////////////////////////////
  Vector(const Vector &obj)
    :mElements(obj.mElements),mComm(obj.mComm),mWorldSize(obj.mWorldSize),mMyRank(obj.mMyRank)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // destructor -- deletes matrix
  //////////////////////////////////////////////////////////////////////////
  ~Vector()
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////
  // additional accessors and accessor prototypes
  //////////////////////////////////////////////////////////////////

  // returns the number of rows
  unsigned int getSize() const { return mElements.size();};


  //////////////////////////////////////////////////////////////////////////
  // v[i] operator -- accessor for elements in vector
  //////////////////////////////////////////////////////////////////////////
  int & operator [](int indx)
  {
    return mElements[indx];
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // v[i] operator -- accessor for elements in vector
  //////////////////////////////////////////////////////////////////////////
  const int & operator [](int indx) const
  {
    return mElements[indx];
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  void resize(int _m)
  {
    mElements.resize(_m);
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // setVal -- function sets all elements to be alpha
  //////////////////////////////////////////////////////////////////////////
  void setScalar(int alpha) 
  {
    for(unsigned int i=0;i<mElements.size();i++)
    {
      mElements[i]=alpha;
    }
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Output vector
  //////////////////////////////////////////////////////////////////////////
  void Print() const
  {
    if(mMyRank==0)
    {
      std::cout << "Vector: " << std::endl;
    }

    for(int rank=0; rank<mWorldSize; rank++)
    {
      if(mMyRank==rank)
      {
        for(unsigned int i=0; i<mElements.size(); i++)
        {
          std::cout << mElements[i] << std::endl;
        }
      }
      MPI_Barrier(mComm);
    }
  }
  //////////////////////////////////////////////////////////////////////////

};
//////////////////////////////////////////////////////////////////////////////

#endif
