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
#include "util.h"
#include "mmio.h"

unsigned int choose2(unsigned int k);

//////////////////////////////////////////////////////////////////////////////
// Enumerate triangles in graph
//////////////////////////////////////////////////////////////////////////////
void Graph::triangleEnumerate()
{
  struct timeval t1, t2;
  double eTime;

  std::cout << "************************************************************"
            << "**********" << std::endl;
  std::cout << "Enumerating triangles ....." << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  std::cout << "--------------------" << std::endl;
  std::cout << "Permuting matrix ...";

  gettimeofday(&t1, NULL);
  mMatrix.permute();
  gettimeofday(&t2, NULL);

  std::cout << " done" <<std::endl;

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);
  std::cout << "TIME - Time to permute  matrix: " << eTime << std::endl;

  std::cout << "--------------------" << std::endl;
  ///////////////////////////////////////////////////////////////////////

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
  // Enumerate triangles
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
  std::cout << "Finished triangle enumeration" << std::endl;
  std::cout << "************************************************************" 
            << "**********" << std::endl;
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Prints triangles in graph
//////////////////////////////////////////////////////////////////////////////
void Graph::printTriangles() const
{
  std::cout << "Triangles: " << std::endl;

  //Iterate through list and output triangles
  std::list<int>::const_iterator iter;
  for (iter=mTriangles.begin(); iter!=mTriangles.end(); iter++)
  {
    std::cout << "(" << *iter+1;
    iter++;
    std::cout << ", " << *iter+1;
    iter++;
    std::cout << ", " << *iter+1 << ")" << std::endl;
  }
}
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//
// Note: Assumes triangles are in order by vertex
//////////////////////////////////////////////////////////////////////////////
void Graph::calculateTriangleDegrees()
{
  std::map<int,int>::iterator vdIter;

  std::map<int,std::map<int,int> >::iterator edIter1;
  std::map<int,int>::iterator edIter2;

  //Iterate through list of triangles
  std::list<int>::const_iterator iter;
  int numEdgesInTriangles=0;
  for (iter=mTriangles.begin(); iter!=mTriangles.end(); iter++)
  {
    int v1 = *iter;
    iter++;
    int v2 = *iter;
    iter++;
    int v3 = *iter;

    /////////////////////////////////
    // Increment triangle degree of v1
    /////////////////////////////////
    vdIter = mVDegrees.find(v1);
    if(vdIter != mVDegrees.end())
    {
      (*vdIter).second++;
    }
    else
    {
      mVDegrees[v1]=1;
    }
    /////////////////////////////////

    /////////////////////////////////
    // Increment triangle degree of v2
    /////////////////////////////////
    vdIter = mVDegrees.find(v2);
    if(vdIter != mVDegrees.end())
    {
      (*vdIter).second++;
    }
    else
    {
      mVDegrees[v2]=1;
    }
    /////////////////////////////////

    /////////////////////////////////
    // Increment triangle degree of v3
    /////////////////////////////////
    vdIter = mVDegrees.find(v3);
    if(vdIter != mVDegrees.end())
    {
      (*vdIter).second++;
    }
    else
    {
      mVDegrees[v3]=1;
    }
    /////////////////////////////////

    /////////////////////////////////
    // Increment triangle degree of edge v1,v2
    /////////////////////////////////
    edIter1 = mEDegrees.find(v1);
    if(edIter1 != mEDegrees.end())
    {
      edIter2 = (*edIter1).second.find(v2);

      if(edIter2 != (*edIter1).second.end())
      {
        (*edIter2).second++;
      }
      else
      {
        (*edIter1).second[v2]=1;
	numEdgesInTriangles++;
      }
    }
    else
    {
      mEDegrees[v1][v2]=1;
      numEdgesInTriangles++;
    }
    /////////////////////////////////

    /////////////////////////////////
    // Increment triangle degree of edge v1,v3
    /////////////////////////////////
    edIter1 = mEDegrees.find(v1);
    if(edIter1 != mEDegrees.end())
    {
      edIter2 = (*edIter1).second.find(v3);

      if(edIter2 != (*edIter1).second.end())
      {
        (*edIter2).second++;
      }
      else
      {
        (*edIter1).second[v3]=1;
        numEdgesInTriangles++;
      }
    }
    else
    {
      mEDegrees[v1][v3]=1;
      numEdgesInTriangles++;
    }
    /////////////////////////////////

    /////////////////////////////////
    // Increment triangle degree of edge v2,v3
    /////////////////////////////////
    edIter1 = mEDegrees.find(v2);
    if(edIter1 != mEDegrees.end())
    {
      edIter2 = (*edIter1).second.find(v3);

      if(edIter2 != (*edIter1).second.end())
      {
        (*edIter2).second++;
      }
      else
      {
        (*edIter1).second[v3]=1;
	numEdgesInTriangles++;
      }
    }
    else
    {
      mEDegrees[v2][v3]=1;
      numEdgesInTriangles++;
    }
    ////////////////////////////////
  }

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Graph::orderTriangles()
{
  int tmpV1, tmpV2, tmpV3;

  //Iterate through list of triangles
  std::list<int>::iterator iter1,iter2,iter3;
  for (iter1=mTriangles.begin(); iter1!=mTriangles.end(); iter1++)
  {
    iter2=iter1; iter2++;
    iter3=iter2; iter3++;

    tmpV1 = *iter1;
    tmpV2 = *iter2;
    tmpV3 = *iter3;

    /////////////////////////////////
    // order triangles
    /////////////////////////////////
    if(tmpV1 < tmpV2)
    {
      if (tmpV1 < tmpV3)
      {
        *iter1 = tmpV1;
	if(tmpV2<tmpV3)
	{
	  *iter2 = tmpV2;
	  *iter3 = tmpV3;
	}
	else
	{
	  *iter2 = tmpV3;
	  *iter3 = tmpV2;
	}
      }
      else
      {
        *iter1 = tmpV3;
        *iter2 = tmpV1;
        *iter3 = tmpV2;
      }

    }
    else
    {
      if (tmpV1 < tmpV3)
      {
        *iter1 = tmpV2;
        *iter2 = tmpV1;
        *iter3 = tmpV3;
      }
      else
      {
        *iter3 = tmpV1;
        if(tmpV2<tmpV3)
        {
          *iter1 = tmpV2;
	  *iter2 = tmpV3;
        }
        else
        {
          *iter1 = tmpV3;
	  *iter2 = tmpV2;
        }
      } 
    }
    /////////////////////////////////

    iter1++;
    iter1++;
  }

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Graph::calculateKCounts()
{
  std::list<int>::const_iterator iter;
  for (iter=mTriangles.begin(); iter!=mTriangles.end(); iter++)
  {
    int v1=*iter;
    iter++;
    int v2=*iter;
    iter++;
    int v3=*iter;

    unsigned int tvMin, teMin;

    findMinTriDegrees(v1,v2,v3,tvMin,teMin);

    int maxK=3;
    for(unsigned int k=3; k<mKCounts.size(); k++)
    {
      if(tvMin >= choose2(k-1) && teMin >= k-2)
      {
        maxK = k;
      }
      else
      {
        break;
      }
    }
    mKCounts[maxK]++;
  }

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Graph::printKCounts()
{
  std::cout << "K-Counts: " << std::endl;
  for(unsigned int i=3; i<mKCounts.size(); i++)
  {
    std::cout << "K[" << i << "] = " <<mKCounts[i] << std::endl;
  }
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Graph::findMinTriDegrees(int v1, int v2, int v3, unsigned int &tvMin, unsigned int &teMin)
{
  /////////////////////////////////////////////////////////////////////////
  // Find tvMin
  /////////////////////////////////////////////////////////////////////////
  tvMin = mVDegrees[v1];
  if(mVDegrees[v2] < (int) tvMin)
  {
    tvMin=mVDegrees[v2];
  }
  if(mVDegrees[v3] < (int) tvMin)
  {
    tvMin=mVDegrees[v3];
  }
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Find teMin
  /////////////////////////////////////////////////////////////////////////
  teMin = mEDegrees[v1][v2];
  if(mEDegrees[v1][v3] < (int) teMin)
  {
    teMin=mEDegrees[v1][v3];
  }
  if(mEDegrees[v2][v3] < (int) teMin)
  {
    teMin=mEDegrees[v2][v3];
  }
  /////////////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
unsigned int choose2(unsigned int k)
{
  if(k==1)
  {
    return 0;
  }
  else if(k>1)
  {
    return k*(k-1)/2;
  }
  return 0;
}
//////////////////////////////////////////////////////////////////////////////

