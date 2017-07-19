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

#include <boost/shared_array.hpp>

//////////////////////////////////////////////////////////////////////////////
// Compressed Sparse Row storage format Matrix
//////////////////////////////////////////////////////////////////////////////
class Vector
{

 private:
  int mNumElements;   //number of elements

  boost::shared_array<int> mElements;                             // elements in vector

  int mBlockSize;

 public:

  //Array of futures
  std::vector<hpx::shared_future<void> > mOps;

  //////////////////////////////////////////////////////////////////////////
  // default constructor -- builds empty matrix
  //////////////////////////////////////////////////////////////////////////
  Vector() 
    :mNumElements(0),mElements(),mBlockSize(1)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // constructor -- allocates memory for CSR sparse matrix
  //////////////////////////////////////////////////////////////////////////
  Vector(int _m, int _bs=1)
    :mNumElements(_m), mElements(new int[_m]()),mBlockSize(_bs)
  {
  };
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Copy constructor -- No data reallocation, just copying of smart pointers
  // TODO: turn mOps into shared_array
  //////////////////////////////////////////////////////////////////////////
  Vector(const Vector &obj)
    :mNumElements(obj.mNumElements),mElements(obj.mElements),mBlockSize(obj.mBlockSize),
    mOps(obj.mOps)
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
  int getSize() const { return mNumElements;};


  //////////////////////////////////////////////////////////////////////////
  // v[i] operator -- accessor for elements in vector
  //////////////////////////////////////////////////////////////////////////
  int operator [](int indx) const
  {
    return mElements[indx];
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  void resize(int _m)
  {
    mNumElements = _m;
    mElements = boost::shared_array<int>(new int[_m]());
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  void setBlockSize(int _bs)
  {
    mBlockSize = _bs;
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // setVal -- function sets element indx to be alpha
  //////////////////////////////////////////////////////////////////////////
  void setVal(int indx,int alpha) 
  {
    mElements[indx]=alpha;
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // setVal -- function sets all elements to be alpha
  //////////////////////////////////////////////////////////////////////////
  void setScalar(int alpha) 
  {
    for(int i=0;i<mNumElements;i++)
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
    std::cout << "Vector: " << std::endl;
    for(int i=0; i<mNumElements; i++)
    {
      std::cout << mElements[i] << std::endl;
    }
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  // Need to optimize this
  //////////////////////////////////////////////////////////////////////////
  void waitOnIndices(const std::set<int> &vecIDs) const
  {
    std::set<int> futIDs;
    std::set<int>::iterator it;

    int currID = 0;
    int sID = 0;
    int eID = sID+mBlockSize;

    // Find set of future indices that need to be waited on

    for (it = vecIDs.begin(); it != vecIDs.end(); ++it)
    {
      if( (*it) < eID )
      {
        futIDs.insert(currID);
      }
      else
      {
        // Do I need to check maxLength of vector?
        // think this should be ok
        while( (*it) >= eID )
	{
          sID = eID;
          eID = sID+mBlockSize;
          currID++;
	}
        futIDs.insert(currID);
      }
    }

    // Wait on all futures that need to be waited on
    for (it = futIDs.begin(); it != futIDs.end(); ++it)
    {
      mOps[(*it)].get();
    }
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  void waitAll() const
  {
    hpx::wait_all(mOps);
  }
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////  
  //////////////////////////////////////////////////////////////////////////
  void clearFutures()
  {
    // perhaps overkill
    std::vector<hpx::shared_future<void> >().swap(mOps);
  }
  //////////////////////////////////////////////////////////////////////////



};
//////////////////////////////////////////////////////////////////////////////

#endif
