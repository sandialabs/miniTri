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
#include "mmio.h"

unsigned int choose2(unsigned int k);


//////////////////////////////////////////////////////////////////////////////
// Enumerate triangles in graph
//////////////////////////////////////////////////////////////////////////////
void Graph::triangleEnumerate()
{
  struct timeval t1, t2;
  double eTime;

  if(mMyRank==0)
  {
    std::cout << "************************************************************"
              << "**********" << std::endl;
    std::cout << "Enumerating triangles ....." << std::endl;
    std::cout << "************************************************************" 
              << "**********" << std::endl;
  }
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  if(mMyRank==0)
  {
    std::cout << "--------------------" << std::endl;
    std::cout << "Permuting matrix ...";
    std::cout << "Permutation code not written!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }

  MPI_Barrier(mComm);
  gettimeofday(&t1, NULL);
  //mMatrix.permute();
  MPI_Barrier(mComm);
  gettimeofday(&t2, NULL);

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);

  if(mMyRank==0)
  {
    std::cout << " done" <<std::endl;

    std::cout << "TIME - Time to permute  matrix: " << eTime << std::endl;

    std::cout << "--------------------" << std::endl;
  }
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  if(mMyRank==0)
  {
    std::cout << "--------------------" << std::endl;
    std::cout << "Creating L,U matrices ...";
  }

  MPI_Barrier(mComm);
  gettimeofday(&t1, NULL);

  CSRMat L(LOWERTRI);
  L.createTriMatrix(mMatrix, LOWERTRI);

  CSRMat U(UPPERTRI);
  U.createTriMatrix(mMatrix, UPPERTRI);

  MPI_Barrier(mComm);
  gettimeofday(&t2, NULL);

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);

  if(mMyRank==0)
  {
    std::cout << " done" <<std::endl;
    std::cout << "TIME - Time to create L, U  matrices: " << eTime << std::endl;

    std::cout << "--------------------" << std::endl;
  }

  //L.print();
  //U.print();

  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // B = L*U
  ///////////////////////////////////////////////////////////////////////
  if(mMyRank==0)
  {
    std::cout << "--------------------" << std::endl;
    std::cout << "B = L*U: " << std::endl;
  }

  CSRMat B(L.getGlobNumRows(),U.getGlobNumCols(),mComm);

  MPI_Barrier(mComm);
  gettimeofday(&t1, NULL);
  B.matmat(L,U);

  MPI_Barrier(mComm);
  gettimeofday(&t2, NULL);

  //B.print();

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);

  if(mMyRank==0)
  {
    std::cout << "TIME - Time to compute B = L*U: " << eTime << std::endl;

    std::cout << "--------------------" << std::endl;
  }
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // B = B .* L
  ///////////////////////////////////////////////////////////////////////
  if(mMyRank==0)
  {
    std::cout << "--------------------" << std::endl;
    std::cout << "B = B .* L: " << std::endl;
  }
  MPI_Barrier(mComm);
  gettimeofday(&t1, NULL);
  B.EWMult(L);
  MPI_Barrier(mComm);
  gettimeofday(&t2, NULL);

  //  B.print();

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);


  if(mMyRank==0)
  {
    std::cout << "TIME - Time to compute B = B .* L: " << eTime << std::endl;
    std::cout << "--------------------" << std::endl;
  }
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // Enumerate triangles
  ///////////////////////////////////////////////////////////////////////
  if(mMyRank==0)
    std::cout << "--------------------" << std::endl;

  MPI_Barrier(mComm);
  gettimeofday(&t1, NULL);
  mTriangles = B.getSumElements();
  gettimeofday(&t2, NULL);

  assert(mTriangles.size() % 3 == 0);
  mLocNumTriangles = mTriangles.size() / 3;

  MPI_Allreduce(&mLocNumTriangles,&mGlobNumTriangles,1,MPI_INT, MPI_SUM, mComm);

  eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);

  if(mMyRank==0)
  {
    std::cout << "TIME - Time to sum up triangles: " << eTime << std::endl;

    std::cout << "--------------------" << std::endl;
  }
  ///////////////////////////////////////////////////////////////////////

  if(mMyRank==0)
  {
    std::cout << "************************************************************"
              << "**********" << std::endl;
    std::cout << "Finished triangle enumeration" << std::endl;
    std::cout << "************************************************************" 
              << "**********" << std::endl;
  }
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
// * Assumes triangles are in order by vertex
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
    /////////////////////////////////

  }

  /////////////////////////////////////////////////////////////
  // Communicate vertex and edge degrees to remote processors that contain
  // triangle with that same vertex
  // -- NOTE: this procedure is inefficient and needs to be rewritten
  /////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////
  // Pack vertex degrees to send
  ////////////////////////////////////////////////////////
  std::vector<int> sendPacked(2*mVDegrees.size());

  unsigned int indx=0;
  for (vdIter=mVDegrees.begin(); vdIter!=mVDegrees.end(); vdIter++)
  {
    sendPacked[indx] = (*vdIter).first;
    indx++;
    sendPacked[indx] = (*vdIter).second;
    indx++;
  }
  ////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////
  // Step by step reduce vertex degrees
  ////////////////////////////////////////////////////////
  MPI_Status status;
  for(int phase=1; phase<mWorldSize; phase++)
  {
    int src = (mMyRank + phase) % mWorldSize;
    int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

    int sendSize = sendPacked.size();
    int recvSize;

    MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0,
                 &recvSize, 1, MPI_INT, src, 0, mComm, &status);


    std::vector<int> recvPacked(recvSize);


    MPI_Sendrecv(&(sendPacked[0]), sendSize, MPI_INT, dst, 0,
                 &(recvPacked[0]), recvSize, MPI_INT, src, 0, mComm, &status);

    // Sum counts if vertex exists in this process's triangles
    for(int i=0;i<recvSize/2;i++)
    {
      vdIter = mVDegrees.find(recvPacked[2*i]);
      if(vdIter != mVDegrees.end())
      {
        (*vdIter).second += recvPacked[2*i+1];
      }
    }

  }
  ////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////
  // Pack edge degrees to send
  ////////////////////////////////////////////////////////
  sendPacked.resize(3*numEdgesInTriangles);

  indx=0;
  for (edIter1=mEDegrees.begin(); edIter1!=mEDegrees.end(); edIter1++)
  {
    for (edIter2=(*edIter1).second.begin(); edIter2!=(*edIter1).second.end(); edIter2++)
    { 
      sendPacked[indx] = (*edIter1).first;
      indx++;
      sendPacked[indx] = (*edIter2).first;
      indx++;
      sendPacked[indx] = (*edIter2).second;
      indx++;
    }
  }
  ////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////
  // Step by step reduce vertex degrees
  ////////////////////////////////////////////////////////
  for(int phase=1; phase<mWorldSize; phase++)
  {
    int src = (mMyRank + phase) % mWorldSize;
    int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

    int sendSize = sendPacked.size();
    int recvSize;

    MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0,
                 &recvSize, 1, MPI_INT, src, 0, mComm, &status);

    std::vector<int> recvPacked(recvSize);

    MPI_Sendrecv(&(sendPacked[0]), sendSize, MPI_INT, dst, 0,
                 &(recvPacked[0]), recvSize, MPI_INT, src, 0, mComm, &status);

    // Sum counts if edge exists in this process's triangles
    for(int i=0;i<recvSize/3;i++)
    {
      edIter1 = mEDegrees.find(recvPacked[3*i]);
      if(edIter1 != mEDegrees.end())
      {
        edIter2 = (*edIter1).second.find(recvPacked[3*i+1]);

        if(edIter2 != (*edIter1).second.end())
        {
          (*edIter2).second += recvPacked[3*i+2];
	}
      }
    }

  }
  ////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////


//   for (vdIter=mVDegrees.begin(); vdIter!=mVDegrees.end(); vdIter++)
//   {
//     std::cout << "V " << (*vdIter).first << ": " << (*vdIter).second << std::endl;
//   }


//   for (edIter1=mEDegrees.begin(); edIter1!=mEDegrees.end(); edIter1++)
//   {
//     for (edIter2=(*edIter1).second.begin(); edIter2!=(*edIter1).second.end(); edIter2++)
//     { 
//       std::cout << "E " << (*edIter1).first << " " << (*edIter2).first << " " << (*edIter2).second << std::endl;
//     }
//   }

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
  std::vector<int> locKCounts(mKCounts.size());


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
    for(unsigned int k=3; k<locKCounts.size(); k++)
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
    locKCounts[maxK]++;
  }


  // Reduce frequency table across processors
  // not sure if src/dst buffs can be same.
  MPI_Allreduce(&(locKCounts[0]),&(mKCounts[0]),locKCounts.size(), MPI_INT, MPI_SUM, mComm);


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
