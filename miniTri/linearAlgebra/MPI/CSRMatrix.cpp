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
// File:      CSRmatrix.cc                                                  //
// Project:   miniTri                                                       //   
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Source file for CSR matrix class.                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <cassert>
#include <cstdlib>

#include "CSRMatrix.hpp"
#include "Vector.hpp"
#include "mmUtil.h"
#include "binFileReader.h"

void printSubmat(const CSRSubmat &submat, int startRow, int locNumRows);

int addNZ(std::map<int,std::list<int> > &nzMap,int col, int elemToAdd);

unsigned int choose2(unsigned int k);

void serialSubmatrixMult(const CSRSubmat & submatA, int numRows,
                         const CSRSubmat & submatB, int startRowB,
                         std::vector<std::map<int,std::list<int> > > & submatCNZs);

void sendRecvSubmat(const CSRSubmat &submatToSend, int dst, int numRowsSend, int startRowSend,
                    MPI_Comm comm,int src, CSRSubmat &remSubmat,int &numRowsRecv,int &recvSubmatStartRow);

//////////////////////////////////////////////////////////////////////////////
// print function -- outputs matrix to file
//                -- accepts optional filename, "CSRmatrix.out" default name
//////////////////////////////////////////////////////////////////////////////
void CSRMat::print() const
{
  std::cout << "Matrix: " << mGlobNumRows << " " << mGlobNumCols << " " << mGlobNNZ  
            << " " << mLocNNZ << std::endl;

  for(int rank=0; rank<mWorldSize; rank++)
  {
    if(mMyRank==rank)
    {
      for(int submatNum=0; submatNum<mWorldSize; submatNum++)
      {
        const CSRSubmat &submat = mSubmat[submatNum];
	printSubmat(submat,mStartRow,mLocNumRows);
      }
    }
    MPI_Barrier(mComm);
  }
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// printSubmat()
//////////////////////////////////////////////////////////////////////////////
void printSubmat(const CSRSubmat &submat, int startRow, int locNumRows)
{
  std::list<int>::const_iterator it;

  for(int rownum=0; rownum<locNumRows; rownum++)
  {
    for(unsigned int nzIdx=0; nzIdx<submat.cols[rownum].size(); nzIdx++)
    {
      std::cout << startRow + rownum << " " << submat.cols[rownum][nzIdx] << " "
                << submat.vals[rownum][nzIdx] << std::endl;
    }
  }

}
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Sums matrix elements
////////////////////////////////////////////////////////////////////////////////
std::list<int> CSRMat::getSumElements() const
{
  std::list<int> matList;

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat &submat = mSubmat[submatNum];
    for(int rownum=0; rownum<mLocNumRows; rownum++)
    {
      int nnz = submat.cols[rownum].size();
      for(int nzIdx=0; nzIdx<nnz; nzIdx++)
      {
        matList.push_back(mStartRow+rownum);
	matList.push_back(submat.vals[rownum][nzIdx]);
	// perhaps should add check for this
	matList.push_back(submat.vals2[rownum][nzIdx]);
      }
    }
  }

  return matList;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// SpMV1 --
//        -- y = this * 1 or y = this' * 1
////////////////////////////////////////////////////////////////////////////////
void CSRMat::SpMV1(bool trans, Vector &y)
{
  ///////////////////////////////////////////////////////////////////////////
  // y = this * 1
  ///////////////////////////////////////////////////////////////////////////
  if(trans==false)
  {
    for (int rowID=0; rowID<mLocNumRows; rowID++)
    {
      y[rowID] = 0;
    } // end loop over rows

    for(unsigned int submatNum=0; submatNum< mSubmat.size(); submatNum++)
    {     
      for (int rowID=0; rowID<mLocNumRows; rowID++)
      {
        y[rowID] += mSubmat[submatNum].cols[rowID].size();
      } // end loop over rows
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // y = this' * 1
  ///////////////////////////////////////////////////////////////////////////
  else
  {

    //////////////////////////////////////////////////////////////////////
    // Purely local computation
    //////////////////////////////////////////////////////////////////////
    const CSRSubmat &submat = mSubmat[mMyRank];
    int colOffset = mSubmatStartCols[mMyRank];

    for (int rowID=0; rowID<mLocNumRows; rowID++)
    {
      std::vector<int>::const_iterator iter;
      for(iter=submat.cols[rowID].begin();iter!=submat.cols[rowID].end();++iter)
      {
        y[(*iter)-colOffset]++;
      }
    } // end loop over rows
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Remote computation
    //////////////////////////////////////////////////////////////////////
    MPI_Status status;

    // Contains remote Y values calculated here
    std::vector<int> sendY;
    std::vector<int> recvY(y.getSize());
    for(int phase=1; phase<mWorldSize; phase++)
    {
      // Will receive data from process src
      int src = (mMyRank + phase) % mWorldSize;
      // performing subcalculation for y elements that belows to process dst
      int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

      //////////////////////////////////////////////////////////////
      // Computation of remote values of y
      //////////////////////////////////////////////////////////////

      // Resize vector and initialize to 0
      sendY.assign(mSubmatNumCols[dst],0);

      const CSRSubmat &submatR = mSubmat[dst];

      int colOffsetR = mSubmatStartCols[dst];
      for (int rowID=0; rowID<mLocNumRows; rowID++)
      {
        std::vector<int>::const_iterator iter;
        for(iter=submatR.cols[rowID].begin();iter!=submatR.cols[rowID].end();++iter)
        {
          sendY[(*iter)-colOffsetR]++;
        }
      } // end loop over rows
      //////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////
      // Communication partial y elements to processor that owns those
      // elements.  Sum into y.
      //////////////////////////////////////////////////////////////
      MPI_Sendrecv(sendY.data(), sendY.size(), MPI_INT, dst, 0,
		   recvY.data(), recvY.size(), MPI_INT, src, 0,
		   mComm, &status);

      for(unsigned int i=0;i<y.getSize();i++)
      {
        y[i] += recvY[i];
      }
      //////////////////////////////////////////////////////////////

    }
    //////////////////////////////////////////////////////////////////////

  }

  return;
}
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// matmat -- level 3 basic linear algebra subroutine
//        -- Z = AB where Z = this
//////////////////////////////////////////////////////////////////////////////
void CSRMat::matmat(const CSRMat &A, const CSRMat &B)
{
  /////////////////////////////////////////////////////////////////////////
  // set dimensions of matrix, build arrays nnzInRow, vals, cols
  /////////////////////////////////////////////////////////////////////////
  int oldGlobNumRows=mGlobNumRows;
  mGlobNumRows = A.getGlobNumRows();
  mGlobNumCols = B.getGlobNumCols();
  mLocNumRows = A.getLocNumRows();
  mStartRow = A.getStartRow();

  // Need to adjust  mSubmatStartCols,mSubmatNumCols if cols change?                                        

  if(oldGlobNumRows!=mGlobNumRows)
  {
    mSubmat.resize(mWorldSize);

    for(unsigned int submatNum=0; submatNum< mSubmat.size(); submatNum++)
    {
      mSubmat[submatNum].cols.resize(mLocNumRows);
      mSubmat[submatNum].vals.resize(mLocNumRows);
      mSubmat[submatNum].vals2.resize(mLocNumRows);

      // This might be unnecessary
      for(int i=0;i<mLocNumRows;i++)
      {
        mSubmat[submatNum].cols[i].resize(0);
	mSubmat[submatNum].vals[i].resize(0);
	mSubmat[submatNum].vals2[i].resize(0);
      }

    }
  }
  else // Allocate space for vals2
  {
    for(unsigned int submatNum=0; submatNum< mSubmat.size(); submatNum++)
    {
      mSubmat[submatNum].vals2.resize(mLocNumRows);
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Structure to hold new nonzeros
  ///////////////////////////////////////////////////////////////////////////
  std::vector<std::vector<std::map<int,std::list<int> > > > CNZs(mWorldSize);
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Completely local computation
  //       Processor 0's perspective: C(1,*) += A(1,1)*B(1,*)
  ///////////////////////////////////////////////////////////////////////////
  const CSRSubmat & submatA = A.getSubMatrix(mMyRank);

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    CNZs[submatNum].resize(mLocNumRows);

    const CSRSubmat & submatB = B.getSubMatrix(submatNum);

    int startRowB = B.getStartRow();

    //////////////////////////////////////////////////////////////////////
    // Updates block: C(myrank,submatNum) += A(myrank,myrank) * B(myrank,submatNum)
    //////////////////////////////////////////////////////////////////////
    serialSubmatrixMult(submatA,A.getLocNumRows(),submatB,startRowB,CNZs[submatNum]);
    //////////////////////////////////////////////////////////////////////
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Communication/Partially remote compuation
  ///////////////////////////////////////////////////////////////////////////
  for(int phase=1; phase<mWorldSize; phase++)
  {
    int src = (mMyRank + phase) % mWorldSize;
    int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

    for(int submatNum=0; submatNum<mWorldSize; submatNum++)
    {
      const CSRSubmat & submatBToSend = B.getSubMatrix(submatNum);

      int remStartRowB=0;
      int remNumRowsB=0;

      CSRSubmat remSubmatB;

      sendRecvSubmat(submatBToSend,dst,B.getLocNumRows(),B.getStartRow(),
		     mComm,src,remSubmatB,remNumRowsB,remStartRowB);

      const CSRSubmat & submatA = A.getSubMatrix(src);

      //////////////////////////////////////////////////////////////////////
      // Updates block: C(myrank,submatNum) += A(myrank,src) * B(src,submatNum)
      //////////////////////////////////////////////////////////////////////
      serialSubmatrixMult(submatA,A.getLocNumRows(),remSubmatB,remStartRowB,CNZs[submatNum]);
      //////////////////////////////////////////////////////////////////////
    }

  }
  ///////////////////////////////////////////////////////////////////////////

  // nonzero stripping should occur here

  //////////////////////////////////////////////////////////////////////
  // Strip out any nonzeros that have only one element
  //   This is an optimization for Triangle Enumeration
  //   Algorithm #2
  //////////////////////////////////////////////////////////////////////
  std::map<int,std::list<int> >::iterator iter;

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    for(int rownum=0;rownum<mLocNumRows;rownum++)
    {
      for (iter=CNZs[submatNum][rownum].begin(); iter!=CNZs[submatNum][rownum].end(); )
      {
        if((*iter).second.size()==1)
        {
          //CNZs[submatNum][rownum].erase(iter++);      // Remove nonzero (C++98 Compliant)
	  iter = CNZs[submatNum][rownum].erase(iter); // Remove nonzero (C++11 Required)
        }
        else
        {
          ++iter;
        }
      }
    }
  }
  //////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Copy resulting nonzeros into matrix object
  ///////////////////////////////////////////////////////////////////////////
  mLocNNZ=0;
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    CSRSubmat & submatC = mSubmat[submatNum];

    for (int rownum=0; rownum<mLocNumRows; rownum++)
    {
      const std::map<int,std::list<int> > &nzMap = CNZs[submatNum][rownum];

      unsigned int nnzInRow = nzMap.size();

      // If there are nonzeros in this row
     if(nnzInRow > 0)
     {
       /////////////////////////////////////////
       // Allocate memory for this row
       /////////////////////////////////////////
       submatC.cols[rownum].resize(nnzInRow);
       submatC.vals[rownum].resize(nnzInRow);
       submatC.vals2[rownum].resize(nnzInRow);

       /////////////////////////////////////////
       // Copy new data into row
       /////////////////////////////////////////
       std::map<int,std::list<int> >::const_iterator iter;
       int nzcnt=0;

       // Iterate through map
       for (iter=nzMap.begin(); iter!=nzMap.end(); iter++)
       {
         submatC.cols[rownum][nzcnt]= (*iter).first;

	 std::list<int>::const_iterator lIter=(*iter).second.begin();
	 submatC.vals[rownum][nzcnt] = *lIter;
	 lIter++;
	 submatC.vals2[rownum][nzcnt]= *lIter;

	 nzcnt++;
       }
       /////////////////////////////////////////
     }

     /////////////////////////////////////////

     mLocNNZ += nnzInRow;

     } // Loop over rows

  } // loop over submats
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  MPI_Allreduce(&mLocNNZ, &mGlobNNZ, 1, MPI_INT, MPI_SUM,mComm);
  ///////////////////////////////////////////////////////////////////////////


}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Multiplies two submatrices together
// Inserts resulting nonzeros into vector of <col,val> value maps
////////////////////////////////////////////////////////////////////////////////
void serialSubmatrixMult(const CSRSubmat & submatA, int numARows,
                         const CSRSubmat & submatB, int startRowB,
                         std::vector<std::map<int,std::list<int> > > & submatCNZs)
{

  ///////////////////////////////////////////////////////////////////////////
  // Compute matrix entries one row at a time
  ///////////////////////////////////////////////////////////////////////////
  for (int rownum=0; rownum<numARows; rownum++)
  {
    int nnzInRowA = submatA.cols[rownum].size();

    for(int nzindxA=0; nzindxA<nnzInRowA; nzindxA++)
    {
      int colA = submatA.cols[rownum][nzindxA];

      int nnzInRowB = submatB.cols[colA-startRowB].size();

      for(int nzindxB=0; nzindxB<nnzInRowB; nzindxB++)
      {
        int colB=submatB.cols[colA-startRowB][nzindxB];

	addNZ(submatCNZs[rownum], colB, colA);
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Send/recv submatrix
////////////////////////////////////////////////////////////////////////////////
void sendRecvSubmat(const CSRSubmat &submatToSend, int dst,
                    int numRowsSend, int startRowSend,
                    MPI_Comm comm, int src,
                    CSRSubmat &remSubmat, int &numRowsRecv, int &startRowRecv)
{
  MPI_Status status;

  /////////////////////////////////////////////////////////////////////////
  // Communicate number of rows in matrix that will be sent
  /////////////////////////////////////////////////////////////////////////
  MPI_Sendrecv(&numRowsSend, 1, MPI_INT, dst, 0,
               &numRowsRecv, 1, MPI_INT, src, 0,
               comm, &status);
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Communicate start row of submatrix to be sent
  /////////////////////////////////////////////////////////////////////////
  MPI_Sendrecv(&startRowSend, 1, MPI_INT, dst, 0,
               &startRowRecv, 1, MPI_INT, src, 0,
               comm, &status);
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Allocate memory for Submatrix struct
  /////////////////////////////////////////////////////////////////////////
  std::vector<int> nnzInRowRecv(numRowsRecv);
  remSubmat.cols.resize(numRowsRecv);
  remSubmat.vals.resize(numRowsRecv);
  std::vector<int> nnzInRowSend(numRowsSend);
  for(unsigned int i=0;i<numRowsSend;i++)
  {
    nnzInRowSend[i] = submatToSend.cols[i].size();
  }
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Sendrecv nnzInRow data
  /////////////////////////////////////////////////////////////////////////
  MPI_Sendrecv(nnzInRowSend.data(), numRowsSend, MPI_INT, dst, 0,
               nnzInRowRecv.data(), numRowsRecv, MPI_INT, src, 0,
               comm, &status);
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // loop over nnzInRow data to allocate cols, vals
  /////////////////////////////////////////////////////////////////////////
  for(int rownum=0; rownum<numRowsRecv; rownum++)
  {
    if(nnzInRowRecv[rownum] > 0)
    {
      remSubmat.cols[rownum].resize(nnzInRowRecv[rownum]);
      remSubmat.vals[rownum].resize(nnzInRowRecv[rownum]);
    }
  }
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Sendrecv col data
  /////////////////////////////////////////////////////////////////////////
  MPI_Request requests[numRowsRecv];
  MPI_Status stats[numRowsRecv];

  for(int rownum=0; rownum<numRowsRecv; rownum++)
  {
    // post irecvs
    // Is this ok when size==0
    if(nnzInRowRecv[rownum] > 0)
    {
      MPI_Irecv((remSubmat.cols[rownum].data()),nnzInRowRecv[rownum], MPI_INT,
		 src, 0, comm, &(requests[rownum]));
    }
    else
    {
      requests[rownum] = MPI_REQUEST_NULL;
    }
  }

  // Is this dangerous
  for(int rownum=0; rownum<numRowsSend; rownum++)
  {
    // post sends
    if(nnzInRowSend[rownum] > 0)
    {
      MPI_Send((void *) submatToSend.cols[rownum].data(), nnzInRowSend[rownum],
	       MPI_INT, dst, 0, comm);
    }
  }

  MPI_Waitall(numRowsRecv,requests,stats);
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Sendrecv val data
  /////////////////////////////////////////////////////////////////////////
  for(int rownum=0; rownum<numRowsRecv; rownum++)
  {
    // post irecvs
    // Is this ok when size==0
    if(nnzInRowRecv[rownum] > 0)
    {
      MPI_Irecv((remSubmat.vals[rownum].data()),nnzInRowRecv[rownum], MPI_INT,
		 src, 0, comm, &(requests[rownum]));
    }
    else
    {
      requests[rownum] = MPI_REQUEST_NULL;
    }
  }

  // Is this dangerous
  for(int rownum=0; rownum<numRowsSend; rownum++)
  {
    // post sends
    if(nnzInRowSend[rownum] > 0)
    {
      MPI_Send((void *)submatToSend.vals[rownum].data(), nnzInRowSend[rownum],
	       MPI_INT, dst, 0, comm);
    }
  }

  MPI_Waitall(numRowsRecv,requests,stats);
  /////////////////////////////////////////////////////////////////////////

}
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CSRMat::readMMMatrix(const char *fname)
{
  //////////////////////////////////////////////////////////////
  // Build edge list from MM file                               
  //////////////////////////////////////////////////////////////
  int numGlobVerts;
  int numLocVerts;
  int startVert;
  std::vector<edge_t> edgeList;


  buildDistEdgeListFromMM(fname, mWorldSize, mMyRank, numGlobVerts, numLocVerts,
                          startVert,edgeList);

  mGlobNumRows = numGlobVerts;
  mLocNumRows = numLocVerts;
  mGlobNumCols = numGlobVerts;
  mStartRow = startVert;
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Sets values for mStartRowOnProc
  ///////////////////////////////////////////////////////////////////////////
  for(int rank=0; rank<mWorldSize; rank++)
  {
    int numrows;
    partitionMatrix(mGlobNumRows,mWorldSize,rank,
        	    numrows,mStartRowOnProc[rank]);
  }
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // Sets values for mSubmatStartCols,mSubmatNumCols                                        
  ///////////////////////////////////////////////////////////////////////////
  for(int rank=0; rank<mWorldSize; rank++)
  {
    partitionMatrix(mGlobNumCols,mWorldSize,rank,
		    mSubmatNumCols[rank],mSubmatStartCols[rank]);
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Allocate memory for matrix                                              
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    mSubmat[submatNum].cols.resize(mLocNumRows);
    mSubmat[submatNum].vals.resize(mLocNumRows);
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // read data from edgelist and insert into temporary sets                  
  ///////////////////////////////////////////////////////////////////////////
  std::vector<std::vector< std::map<int,int> > >rowSets(mWorldSize);
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    rowSets[submatNum].resize(mLocNumRows);
  }

  // read in matrix entries from file
  int tmpval=1.0; //Is this necessary

  unsigned int locNNZ = edgeList.size();

  for (unsigned int i=0; i<locNNZ; i++)
  {
    int tmprow = edgeList[i].v0;
    int tmpcol = edgeList[i].v1;

    int subMat = whichSubMatrix(tmpcol-1);

    rowSets[subMat][tmprow-1-mStartRow][tmpcol-1] = tmpval;
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  //clear edgelist                                                           
  ///////////////////////////////////////////////////////////////////////////
  std::vector<edge_t>().swap(edgeList);
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // copy data from temporary sets to matrix data structures                 
  ///////////////////////////////////////////////////////////////////////////
  std::map<int,int>::iterator iter;

  mLocNNZ = 0;

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    std::vector< std::map<int,int> > &rowSetsRef = rowSets[submatNum];
    CSRSubmat &submat = mSubmat[submatNum];

    for(int rownum=0; rownum<mLocNumRows; rownum++)
    {
      int nnzIndx=0;
      int nnzToAdd = rowSetsRef[rownum].size();
      mLocNNZ += nnzToAdd;

      if(nnzToAdd>0)
      {
        submat.cols[rownum].resize(nnzToAdd);
	submat.vals[rownum].resize(nnzToAdd);

	for (iter=rowSetsRef[rownum].begin();iter!=rowSetsRef[rownum].end();iter++)
	{
	  submat.cols[rownum][nnzIndx] = (*iter).first;
	  submat.vals[rownum][nnzIndx] = (*iter).second;
	  nnzIndx++;
	}
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Sum mLocNNZ over all processes to obtain mGlobNNZ                       
  ///////////////////////////////////////////////////////////////////////////
  MPI_Allreduce(&mLocNNZ,&mGlobNNZ,1,MPI_INT, MPI_SUM, mComm);
  ///////////////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CSRMat::createIncidenceMatrix(const CSRMat &adjMatrix, std::map<int,std::map<int,int> > & eIndices)
{

   assert(mSubmat.size()==mWorldSize);

   mGlobNumRows = adjMatrix.getGlobNumRows();
   mLocNumRows  = adjMatrix.getLocNumRows();
   mStartRow = adjMatrix.getStartRow();

  ///////////////////////////////////////////////////////////////////////////
  // Given local nonzeros in adjacency matrix, generate nonzeros for incidence matrix
  ///////////////////////////////////////////////////////////////////////////

  // using tmpVector since we have to update edge IDs based on edges built on
  // remote processors.  Updating set values is a fairly expensive process, so
  // instead we temporarily use vectors and copy into sets once we have the
  // the edge offsets
  std::vector<std::vector<int> > tmpColsInRow(mLocNumRows);

  std::vector<std::vector<int> > remNZs(mWorldSize);

  int eCnt=0;

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat &submat = adjMatrix.getSubMatrix(submatNum);

    for(int rownum=0; rownum<mLocNumRows; rownum++)
    {
      for(unsigned int nzIdx=0; nzIdx<submat.cols[rownum].size(); nzIdx++)
      {
        int colnum = submat.cols[rownum][nzIdx];
        if(rownum+mStartRow < colnum)      
	{
	  tmpColsInRow[rownum].push_back(eCnt);

          // These can be remote and may need to be communicated, submats can probably be used here
          if(submatNum==mMyRank)
	  {
	    tmpColsInRow[colnum-mStartRow].push_back(eCnt); 
          }
          else
	  {
            remNZs[submatNum].push_back(colnum);
            remNZs[submatNum].push_back(eCnt);
	  }

	  eIndices[mStartRow+rownum][colnum]=eCnt;    // Need to handle globally?
          eCnt++;
	}
      }
    }

  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Determine edge number offset using prefix_sum
  // Determine total number of edges/columns using all_reduce
  ///////////////////////////////////////////////////////////////////////////
  int eOffset;
  MPI_Exscan(&eCnt, &eOffset, 1, MPI_INT, MPI_SUM, mComm);

  if(mMyRank==0)
  {
    eOffset=0;
  }

  MPI_Allreduce(&eCnt, &mGlobNumCols, 1, MPI_INT, MPI_SUM, mComm);
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Adjust all edge numbers based on edge offset
  ///////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Adjust vertex pair to edge map
  //////////////////////////////////////////////////////////////////////
  std::map<int,int>::iterator mapIter;

  for(int rownum=mStartRow; rownum<mStartRow+mLocNumRows; rownum++)
  {
    for (mapIter=eIndices[rownum].begin();mapIter!=eIndices[rownum].end();mapIter++)
    {
      (*mapIter).second = (*mapIter).second+eOffset;
    }
  }
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Adjust edge column ID for local nonzeros
  //////////////////////////////////////////////////////////////////////
  std::vector<std::set<int> > colsInRow(mLocNumRows);

  std::vector<int>::const_iterator iter;
  for(int rownum=0; rownum<mLocNumRows; rownum++)
  {
    for (iter=tmpColsInRow[rownum].begin();iter!=tmpColsInRow[rownum].end();iter++)
    {
      colsInRow[rownum].insert( *iter + eOffset);   
    }

    // free memory of vector
    std::vector<int>().swap(tmpColsInRow[rownum]);
  }

  // free memory of vector
  std::vector<std::vector<int> >().swap(tmpColsInRow);
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Adjust edge column ID for remote nonzeros
  //////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    unsigned int nnz = remNZs[submatNum].size() / 2;
 
    for(unsigned int nzIdx=0; nzIdx<nnz; nzIdx++)
    {
      //Update edge column ID
      remNZs[submatNum][2*nzIdx+1] += eOffset;
    }

  }
  //////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Communicate remote nonzeros to process that owns them
  ///////////////////////////////////////////////////////////////////////////
  MPI_Status status;
  std::vector<int> recvBuffer;
  for(int phase=1; phase<mWorldSize; phase++)
  {
    int src = (mMyRank + phase) % mWorldSize;
    int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

    int sendSize = remNZs[dst].size();
    int recvSize;

    MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0, &recvSize, 1, MPI_INT, src, 0,
		 mComm, &status);

    // Need to send data
    if(sendSize!=0)
    {
      // Need to receive data
      if(recvSize!=0)
      {
        recvBuffer.resize(recvSize);
        MPI_Sendrecv(remNZs[dst].data(), sendSize, MPI_INT, dst, 0, 
		     recvBuffer.data(), recvSize, MPI_INT, src, 0,
		     mComm, &status);

        unsigned int nnz=recvSize/2;

        for(unsigned int nzIdx=0; nzIdx<nnz; nzIdx++)
	{
          int rIdx= recvBuffer[2*nzIdx];
          int cIdx= recvBuffer[2*nzIdx+1];
          colsInRow[rIdx-mStartRow].insert(cIdx);
	}
      }
      else
      {
        MPI_Send(remNZs[dst].data(), sendSize, MPI_INT, dst, 0, mComm);
      }

    }
    // Need to receive data
    else if(recvSize!=0)
    {
      recvBuffer.resize(recvSize);
      MPI_Recv(recvBuffer.data(), recvSize, MPI_INT, src, 0, mComm, &status);

      unsigned int nnz=recvSize/2;

      for(unsigned int nzIdx=0; nzIdx<nnz; nzIdx++)
      {
        int rIdx= recvBuffer[2*nzIdx];
        int cIdx= recvBuffer[2*nzIdx+1];
        colsInRow[rIdx-mStartRow].insert(cIdx);
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Sets values for mStartRowOnProc
  ///////////////////////////////////////////////////////////////////////////
  for(int rank=0; rank<mWorldSize; rank++)
  {
    int numrows;
    partitionMatrix(mGlobNumRows,mWorldSize,rank,
         	    numrows,mStartRowOnProc[rank]);
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Sets values for mSubmatStartCols,mSubmatNumCols                                        
  ///////////////////////////////////////////////////////////////////////////
  for(int rank=0; rank<mWorldSize; rank++)
  {
      partitionMatrix(mGlobNumCols,mWorldSize,rank,
	  	      mSubmatNumCols[rank],mSubmatStartCols[rank]);
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Resize memory for matrix 
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    mSubmat[submatNum].cols.resize(mLocNumRows);
    mSubmat[submatNum].vals.resize(mLocNumRows);
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Copy data into matrix data structures
  //     -- Can probably free tmp data structures throughout
  ///////////////////////////////////////////////////////////////////////////
  for(int rownum=0; rownum<mLocNumRows; rownum++)
  {
    int nnzToAdd = colsInRow[rownum].size();

    mLocNNZ += nnzToAdd;

    std::set<int>::const_iterator iter;
    for (iter=colsInRow[rownum].begin();iter!=colsInRow[rownum].end();iter++)
    {
      int col = (*iter);

      int submatNum = whichSubMatrix(col);

      mSubmat[submatNum].cols[rownum].push_back(col);
      mSubmat[submatNum].vals[rownum].push_back(1);
    }
  }
  ///////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////
  // Sum mLocNNZ over all processes to obtain mGlobNNZ                       
  ///////////////////////////////////////////////////////////////////////////
  MPI_Allreduce(&mLocNNZ,&mGlobNNZ,1,MPI_INT, MPI_SUM, mComm);
  ///////////////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// addNZ -- For a given row, add a column for a nonzero into a sorted list
//////////////////////////////////////////////////////////////////////////////
int addNZ(std::map<int,std::list<int> > &nzMap,int col, int elemToAdd)
{
  std::map<int,std::list<int> >::iterator it;

  it = nzMap.find(col);

  //////////////////////////////////////
  //If columns match, no additional nz, add element to end of list
  //////////////////////////////////////      
  if(it != nzMap.end())
  {
    (*it).second.push_back(elemToAdd);
    return 0;
  }

  std::list<int> newList;
  newList.push_back(elemToAdd);
  nzMap.insert(std::pair<int,std::list<int> >(col, newList));
  return 1;
}
//////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Compute K counts
//
////////////////////////////////////////////////////////////////////////////////
void CSRMat::computeKCounts(const Vector &vTriDegrees,const Vector &eTriDegrees,
                            std::map<int,std::map<int,int> > & edgeInds,
                            std::vector<int> &kCounts)
{
  std::vector<int> locKCounts(kCounts.size());

  // Map of vertices in local triangles to vertex triangle degrees
  std::map<int,int> tvMap;

  // Map of pairs of vertices in local triangles to vertex edge degrees
  std::map<int,std::map<int,int> > teMap;


  // Vertex triangle degrees for vertices in local triangles that are owned
  // by remote processors
  std::vector<std::vector<int> > remV(mWorldSize);

  ///////////////////////////////////////////////////////////////////////////
  // 1. If vertices are local, insert triangle degree info into tvMap
  // 2. Otherwise, store in vertices in remV array
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat &submat = mSubmat[submatNum];

    for (int rownum=0; rownum<mLocNumRows; rownum++)
    {
      unsigned int nnz=submat.cols[rownum].size();

      for(int nzIdx=0; nzIdx<nnz; nzIdx++)
      {

        int v1 = mStartRow + rownum;
        int v2 = submat.vals[rownum][nzIdx];
        int v3 = submat.vals2[rownum][nzIdx];

        // Removes extra work from redundant triangles
        if(v1>v2 && v1>v3)
        {
          // local vertex
          tvMap[v1] = vTriDegrees[rownum];

          int p2 = rowOnWhichProc(v2);
          if(mMyRank==p2) // local vertex
	  {
            tvMap[v2] = vTriDegrees[v2-mStartRow];
	  }
          else // remote vertex
	  {
            remV[p2].push_back(v2);
	  }

          int p3 = rowOnWhichProc(v3);
          if(mMyRank==p3) // local vertex
	  {
            tvMap[v3] = vTriDegrees[v3-mStartRow];
	  }
          else // remote vertex
	  {
            remV[p3].push_back(v3);
	  }
	}
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Communicate Tv to remote processes that need it
  ///////////////////////////////////////////////////////////////////////////
  std::vector<int> recvV;
  std::vector<int> recvTv;
  MPI_Status status;
  for(int phase=1; phase<mWorldSize; phase++)
  {
    int src = (mMyRank + phase) % mWorldSize;
    int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

    int sendSize = remV[dst].size();
    int recvSize;

    MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0, &recvSize, 1, MPI_INT, src, 0,
		 mComm, &status);

    // Need to send V data
    if(sendSize!=0)
    {
      // Need to receive V data
      if(recvSize!=0)
      {
        ////////////////////////////////////////////////////////////
        // Sendrecv vertices that belong to remote process
        ////////////////////////////////////////////////////////////
        recvV.resize(recvSize);
	MPI_Sendrecv(remV[dst].data(), sendSize, MPI_INT, dst, 0,
		     recvV.data(), recvSize, MPI_INT, src, 0,
		     mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Fill buffer with Tv corresponding to vertices
        ////////////////////////////////////////////////////////////
        for(int indx=0; indx<recvSize; indx++)
	{
          recvV[indx] = vTriDegrees[recvV[indx]-mStartRow];
        } 
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Sendrecv corresponding Tv values
        ////////////////////////////////////////////////////////////
        recvTv.resize(sendSize);
       
	MPI_Sendrecv(recvV.data(), recvSize, MPI_INT, src, 0,
		     recvTv.data(), sendSize, MPI_INT, dst, 0,
		     mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Insert v/Tv pairs into map
        ////////////////////////////////////////////////////////////
        for(int indx=0; indx<sendSize; indx++)
	{
          tvMap[remV[dst][indx]] = recvTv[indx];
        } 
        //////////////////////////////////////////////////////////// 
      }
      else
      {
        ////////////////////////////////////////////////////////////
        // Send vertices that belong to remote process
        ////////////////////////////////////////////////////////////
        MPI_Send(remV[dst].data(), sendSize, MPI_INT, dst, 0, mComm);

        ////////////////////////////////////////////////////////////
        // 1. Recv corresponding Tv values
        // 2. Insert v/Tv pairs into map
        ////////////////////////////////////////////////////////////
        recvTv.resize(sendSize);
	MPI_Recv(recvTv.data(), sendSize, MPI_INT, dst, 0,
		 mComm, &status);

        for(int indx=0; indx<sendSize; indx++)
	{
          tvMap[remV[dst][indx]] = recvTv[indx];
        } 
        ////////////////////////////////////////////////////////////
      }
    }
    // Need to receive data
    else if(recvSize!=0)
    {
      ////////////////////////////////////////////////////////////
      // Recv vertices that belong to this process
      ////////////////////////////////////////////////////////////
      recvV.resize(recvSize);
      MPI_Recv(recvV.data(), recvSize, MPI_INT, src, 0,
	       mComm, &status);
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // 1. Fill buffer with Tv corresponding to vertices
      // 2. Send corresponding Tv values to proc that needs them
      ////////////////////////////////////////////////////////////
      for(int indx=0; indx<recvSize; indx++)
      {
        recvV[indx] = vTriDegrees[recvV[indx]-mStartRow];
      } 
      MPI_Send(recvV.data(), recvSize, MPI_INT, src, 0, mComm);
      ////////////////////////////////////////////////////////////
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  //CLEAN UP DATA


  ///////////////////////////////////////////////////////////////////////////
  // If vertex pair -> edgeID lookup is local 
  //     If edge is local, insert triangle degree info into teMap
  //     Else store edgeID in remEIDs and edge info in remEdges
  // Else store vertex pair in eIDtoFind 
  ///////////////////////////////////////////////////////////////////////////

  // Vertex pairs that are owned remotely
  // Communication needed to find corresponding edge ids
  std::vector<std::vector<int> > eIDtoFind(mWorldSize);

  std::vector<std::vector<int> > remEIDs(mWorldSize);
  std::vector<std::vector<edge_t> > remEdges(mWorldSize);

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat &submat = mSubmat[submatNum];

    for (int rownum=0; rownum<mLocNumRows; rownum++)
    {
      unsigned int nnz=submat.cols[rownum].size();

      for(int nzIdx=0; nzIdx<nnz; nzIdx++)
      {
        int v1 = mStartRow + rownum;
        int v2 = submat.vals[rownum][nzIdx];
        int v3 = submat.vals2[rownum][nzIdx];

        // Removes redundant triangles 
        // (perhaps should rework so v1 is least, would make v1,v2 -> edge map more local)
        if(v1>v2 && v1>v3)
        {

          /////////////////////////////////////////////////////////////////
          // Vertex order: v2, v3, v1
          /////////////////////////////////////////////////////////////////        
 	  if(v2<v3)
	  {

	    int p2 = rowOnWhichProc(v2);
	    if(mMyRank==p2) // local vertex
	    {
              // edge v2,v3
              int eID1 = edgeInds.find(v2)->second.find(v3)->second;
              int submat1 = whichSubMatrix(eID1);
              if(submat1 == mMyRank)
	      {
                teMap[v2][v3] = eTriDegrees[eID1-mSubmatStartCols[submat1]];
              }
              else
	      {
                remEIDs[submat1].push_back(eID1);

                edge_t e1;
                e1.v0=v2;
                e1.v1=v3;
                remEdges[submat1].push_back(e1);
	      }

              // edge v2,v1
 	      int eID2 = edgeInds.find(v2)->second.find(v1)->second;
              int submat2 = whichSubMatrix(eID2);
              if(submat2 == mMyRank)
	      {
                teMap[v2][v1] = eTriDegrees[eID2-mSubmatStartCols[submat2]];
	      }
              else
	      {
                remEIDs[submat2].push_back(eID2);

                edge_t e2;
                e2.v0=v2;
                e2.v1=v1;
                remEdges[submat2].push_back(e2);
	      }

	    }
	    else // remote vertex
	    {
              // edge v2,v3
              eIDtoFind[p2].push_back(v2);
	      eIDtoFind[p2].push_back(v3);

              // edge v2,v1
              eIDtoFind[p2].push_back(v2);
	      eIDtoFind[p2].push_back(v1);
	    }

            int p3 = rowOnWhichProc(v3);
            if(p3==mMyRank)
	    {
              // edge v3,v1
              int eID3 = edgeInds.find(v3)->second.find(v1)->second;

              int submat3 = whichSubMatrix(eID3);
              if(submat3 == mMyRank)
	      {
                teMap[v3][v1] = eTriDegrees[eID3-mSubmatStartCols[submat3]];
	      }
              else
	      {
                remEIDs[submat3].push_back(eID3);

                edge_t e3;
                e3.v0=v3;
                e3.v1=v1;
                remEdges[submat3].push_back(e3);
	      }

	    }
            else // remote vertex
	    {
              // edge v3,v1
              eIDtoFind[p3].push_back(v3);
	      eIDtoFind[p3].push_back(v1);
            }

	  }
          /////////////////////////////////////////////////////////////////
          // Vertex order: v3, v2, v1
          /////////////////////////////////////////////////////////////////
          else
	  {
	    int p3 = rowOnWhichProc(v3);
	    if(mMyRank==p3) // local vertex
	    {
              // edge v3,v2
              int eID1 = edgeInds.find(v3)->second.find(v2)->second;
              int submat1 = whichSubMatrix(eID1);
              if(submat1 == mMyRank)
	      {
                teMap[v3][v2] = eTriDegrees[eID1-mSubmatStartCols[submat1]];
              }
              else
	      {
                remEIDs[submat1].push_back(eID1); 

                edge_t e1;
                e1.v0=v3;
                e1.v1=v2;
                remEdges[submat1].push_back(e1);
	      }

              // edge v3,v1
 	      int eID2 = edgeInds.find(v3)->second.find(v1)->second;
              int submat2 = whichSubMatrix(eID2);
              if(submat2 == mMyRank)
	      {
                teMap[v3][v1] = eTriDegrees[eID2-mSubmatStartCols[submat2]];
              }
              else
	      {
                remEIDs[submat2].push_back(eID2); 

                edge_t e2;
                e2.v0=v3;
                e2.v1=v1;
                remEdges[submat2].push_back(e2);
	      }
	    }
	    else // remote vertex
	    {
              // edge v3,v2
              eIDtoFind[p3].push_back(v3);
	      eIDtoFind[p3].push_back(v2);

              // edge v3,v1
              eIDtoFind[p3].push_back(v3);
	      eIDtoFind[p3].push_back(v1);
	    }

            int p2 = rowOnWhichProc(v2);
            if(p2==mMyRank)
	    {
              // edge v2,v1
              int eID3 = edgeInds.find(v2)->second.find(v1)->second;

              int submat3 = whichSubMatrix(eID3);
              if(submat3 == mMyRank)
	      {
                teMap[v2][v1] = eTriDegrees[eID3-mSubmatStartCols[submat3]];
              }
              else
	      {
                remEIDs[submat3].push_back(eID3); 

                edge_t e3;
                e3.v0=v2;
                e3.v1=v1;
                remEdges[submat3].push_back(e3);
	      }
	    }
            else
	    {
              // edge v2,v1
              eIDtoFind[p2].push_back(v2);
	      eIDtoFind[p2].push_back(v1);
            }

	  }
          /////////////////////////////////////////////////////////////////

	}
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Communicate edge IDs to remote processes that need it
  ///////////////////////////////////////////////////////////////////////////
  std::vector<int> recvVPairs;

  std::vector<int> sendEIDs;
  std::vector<int> recvEIDs;
  for(int phase=1; phase<mWorldSize; phase++)
  {
    int src = (mMyRank + phase) % mWorldSize;
    int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

    int sendSize = eIDtoFind[dst].size();
    int recvSize;

    // Determine amount of data that needs to be received
    MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0, &recvSize, 1, MPI_INT, src, 0,
		 mComm, &status);

    // Need to send eID data
    if(sendSize!=0)
    {
      // Need to receive eID data
      if(recvSize!=0)
      {
        ////////////////////////////////////////////////////////////
        // Sendrecv eIDs that belong to remote process
        ////////////////////////////////////////////////////////////
        recvVPairs.resize(recvSize);
	MPI_Sendrecv(eIDtoFind[dst].data(), sendSize, MPI_INT, dst, 0,
		     recvVPairs.data(), recvSize, MPI_INT, src, 0,
		     mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Fill buffer with eIDs corresponding to vertex pairs
        ////////////////////////////////////////////////////////////
        int sendSize2 = recvSize/2;
        sendEIDs.resize(sendSize2);
        for(int indx=0; indx<sendSize2; indx++)
	{
          int v1 = recvVPairs[2*indx];
          int v2 = recvVPairs[2*indx+1];
          sendEIDs[indx] = edgeInds.find(v1)->second.find(v2)->second;
        } 
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Sendrecv corresponding EIDs
        ////////////////////////////////////////////////////////////
        int recvSize2 = sendSize/2;
        recvEIDs.resize(recvSize2);
       
	MPI_Sendrecv(sendEIDs.data(), sendSize2, MPI_INT, src, 0,
		     recvEIDs.data(), recvSize2, MPI_INT, dst, 0,
		     mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // If eID is local, insert [v1, v2], Te into map
        // Otherwise, store info in remEIDs and remEdges
        ////////////////////////////////////////////////////////////
        for(int indx=0; indx<recvSize2; indx++)
	{
	  int submat = whichSubMatrix(recvEIDs[indx]);
	  if(submat == mMyRank)
	  {
            int v1 = eIDtoFind[dst][2*indx];
            int v2 = eIDtoFind[dst][2*indx+1];

	    teMap[v1][v2] = eTriDegrees[recvEIDs[indx]-mSubmatStartCols[submat]];
	  }
	  else
	  {
	    remEIDs[submat].push_back(recvEIDs[indx]);

	    edge_t e;
	    e.v0=eIDtoFind[dst][2*indx];
	    e.v1=eIDtoFind[dst][2*indx+1];
	    remEdges[submat].push_back(e);
	  }
        }
        //////////////////////////////////////////////////////////// 
      }
      else
      {
        ////////////////////////////////////////////////////////////
        // Send eIDs that belong to remote process
        ////////////////////////////////////////////////////////////
	MPI_Send(eIDtoFind[dst].data(), sendSize, MPI_INT, dst, 0,
		 mComm);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Recv corresponding EIDs
        ////////////////////////////////////////////////////////////
        int recvSize2 = sendSize/2;
        recvEIDs.resize(recvSize2);
       
	MPI_Recv(recvEIDs.data(), recvSize2, MPI_INT, dst, 0,
		 mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // If eID is local, insert [v1, v2], Te into map
        // Otherwise, store info in remEIDs and remEdges
        ////////////////////////////////////////////////////////////
        for(int indx=0; indx<recvSize2; indx++)
	{
	  int submat = whichSubMatrix(recvEIDs[indx]);
	  if(submat == mMyRank)
	  {
            int v1 = eIDtoFind[dst][2*indx];
            int v2 = eIDtoFind[dst][2*indx+1];

	    teMap[v1][v2] = eTriDegrees[recvEIDs[indx]-mSubmatStartCols[submat]];
	  }
	  else
	  {
	    remEIDs[submat].push_back(recvEIDs[indx]);

	    edge_t e;
	    e.v0=eIDtoFind[dst][2*indx];
	    e.v1=eIDtoFind[dst][2*indx+1];
	    remEdges[submat].push_back(e);
	  }
        }
        //////////////////////////////////////////////////////////// 
      }
    }
    // Need to receive data
    else if(recvSize!=0)
    {
      ////////////////////////////////////////////////////////////
      // Recv eIDs that belong to remote process
      ////////////////////////////////////////////////////////////
      recvVPairs.resize(recvSize);
      MPI_Recv(recvVPairs.data(), recvSize, MPI_INT, src, 0,
	       mComm, &status);
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // Fill buffer with EIDs corresponding to vertex pairs
      ////////////////////////////////////////////////////////////
      int sendSize2 = recvSize/2;
      sendEIDs.resize(sendSize2);
      for(int indx=0; indx<sendSize2; indx++)
      {
        int v1 = recvVPairs[2*indx];
        int v2 = recvVPairs[2*indx+1];
        sendEIDs[indx] = edgeInds.find(v1)->second.find(v2)->second;
      } 
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // Send corresponding EIDs
      ////////////////////////////////////////////////////////////
      MPI_Send(sendEIDs.data(), sendSize2, MPI_INT, src, 0,
	       mComm);
      ////////////////////////////////////////////////////////////

    }
  }
  ///////////////////////////////////////////////////////////////////////////

//CLEAN UP DATA -- remove eIDtoFind, sendEIDs

  ///////////////////////////////////////////////////////////////////////////
  // Communicate edge triangle degrees to remote processes that need it
  ///////////////////////////////////////////////////////////////////////////
  std::vector<int> sendEDegrees;
  std::vector<int> recvEDegrees;
  for(int phase=1; phase<mWorldSize; phase++)
  {
    int src = (mMyRank + phase) % mWorldSize;
    int dst = (mMyRank + mWorldSize-phase) % mWorldSize;

    int sendSize = remEIDs[dst].size();
    int recvSize;

    // Determine amount of data that needs to be received
    MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0, &recvSize, 1, MPI_INT, src, 0,
		 mComm, &status);

    // Need to send eIDs
    if(sendSize!=0)
    {
      // Need to receive eIDs
      if(recvSize!=0)
      {
        ////////////////////////////////////////////////////////////
        // Sendrecv eIDs that belong to remote process
        ////////////////////////////////////////////////////////////
        recvEIDs.resize(recvSize);
	MPI_Sendrecv(remEIDs[dst].data(), sendSize, MPI_INT, dst, 0,
		     recvEIDs.data(), recvSize, MPI_INT, src, 0,
		     mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Fill buffer with E triangle degrees corresponding to EIDs
        ////////////////////////////////////////////////////////////
	std::vector<int> &sendEDegrees = recvEIDs; 
        int sendSize2 = recvSize;
        for(int indx=0; indx<sendSize2; indx++)
	{
	  sendEDegrees[indx] = eTriDegrees[recvEIDs[indx]-mSubmatStartCols[mMyRank]];
        } 
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Sendrecv corresponding E triangle degrees
        ////////////////////////////////////////////////////////////
        int recvSize2 = sendSize;
        recvEDegrees.resize(recvSize2);
       
	MPI_Sendrecv(sendEDegrees.data(), sendSize2, MPI_INT, src, 0,
		     recvEDegrees.data(), recvSize2, MPI_INT, dst, 0,
		     mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Insert [v1, v2], Te into map
        ////////////////////////////////////////////////////////////
        for(int indx=0; indx<recvSize2; indx++)
	{
          int v1 = remEdges[dst][indx].v0;
          int v2 = remEdges[dst][indx].v1;

	  teMap[v1][v2] = recvEDegrees[indx];
        }
        //////////////////////////////////////////////////////////// 
      }
      else
      {
        ////////////////////////////////////////////////////////////
        // Send eIDs that belong to remote process
        ////////////////////////////////////////////////////////////
	MPI_Send(remEIDs[dst].data(), sendSize, MPI_INT, dst, 0,
		 mComm);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Recv corresponding edge triangle degrees
        ////////////////////////////////////////////////////////////
        int recvSize2 = sendSize;
        recvEDegrees.resize(recvSize2);
       
	MPI_Recv(recvEDegrees.data(), recvSize2, MPI_INT, dst, 0,
		 mComm, &status);
        ////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////
        // Insert [v1, v2], Te into map
        // Otherwise, store info in remEID and remEdges
        ////////////////////////////////////////////////////////////
        for(int indx=0; indx<recvSize2; indx++)
	{
          int v1 = remEdges[dst][indx].v0;
          int v2 = remEdges[dst][indx].v1;

	  teMap[v1][v2] = recvEDegrees[indx];
        }
        //////////////////////////////////////////////////////////// 
      }
    }
    // Need to receive data
    else if(recvSize!=0)
    {
      ////////////////////////////////////////////////////////////
      // Recv eIDs that belong to remote process
      ////////////////////////////////////////////////////////////
      recvEIDs.resize(recvSize);
      MPI_Recv(recvEIDs.data(), recvSize, MPI_INT, src, 0, mComm, &status);
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // Fill buffer with E triangle degrees corresponding to EIDs
      ////////////////////////////////////////////////////////////
      std::vector<int> &sendEDegrees = recvEIDs; 
      int sendSize2 = recvSize;
      for(int indx=0; indx<sendSize2; indx++)
      {
	sendEDegrees[indx] = eTriDegrees[recvEIDs[indx]-mSubmatStartCols[mMyRank]];
      } 
      ////////////////////////////////////////////////////////////

      ////////////////////////////////////////////////////////////
      // Send corresponding E triangle degrees
      ////////////////////////////////////////////////////////////
      MPI_Send(sendEDegrees.data(), sendSize2, MPI_INT, src, 0, mComm);
      ////////////////////////////////////////////////////////////
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Loop over all triangles on process, compute kcount
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat &submat = mSubmat[submatNum];

    for (int rownum=0; rownum<mLocNumRows; rownum++)
    {
      unsigned int nnz=submat.cols[rownum].size();

      for(int nzIdx=0; nzIdx<nnz; nzIdx++)
      {
        int v1 = mStartRow + rownum;
        int v2 = submat.vals[rownum][nzIdx];
        int v3 = submat.vals2[rownum][nzIdx];

        // Removes extra work from redundant triangles
        if(v1>v2 && v1>v3)
        {

          ////////////////////////////////////////////////////////////////////
	  // Find tvMin
	  ////////////////////////////////////////////////////////////////////
	  int vDegree1 = tvMap.find(v1)->second;
	  int vDegree2 = tvMap.find(v2)->second;
	  int vDegree3 = tvMap.find(v3)->second;

	  // BUG FIX
// 	  unsigned int tvMin = std::min(std::min(vTriDegrees[v1],vTriDegrees[v2]),
// 					vTriDegrees[v3]);
	  unsigned int tvMin = std::min(std::min(vDegree1,vDegree2),
					vDegree3);
	  ////////////////////////////////////////////////////////////////////


	  /////////////////////////////////////////////////////////////////////////
	  // Find teMin                                                            
	  /////////////////////////////////////////////////////////////////////////

	  // I believe that v2<v3 by construction
	  int eDegree1,eDegree2,eDegree3;
	  if(v2<v3)
	  {
	    eDegree1 = teMap.find(v2)->second.find(v3)->second;
	    eDegree2 = teMap.find(v2)->second.find(v1)->second;
	    eDegree3 = teMap.find(v3)->second.find(v1)->second;
	  }
	  else
	  {
	    eDegree1 = teMap.find(v3)->second.find(v2)->second;
	    eDegree2 = teMap.find(v3)->second.find(v1)->second;
	    eDegree3 = teMap.find(v2)->second.find(v1)->second;
	  }

 	  unsigned int teMin = std::min(std::min(eDegree1,eDegree2),eDegree3);
	  /////////////////////////////////////////////////////////////////////////

	  /////////////////////////////////////////////////////////////////////////
	  // Determine k count for triangle                                        
	  /////////////////////////////////////////////////////////////////////////
	  unsigned int maxK=3;
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
	  /////////////////////////////////////////////////////////////////////////
	}
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  MPI_Allreduce(locKCounts.data(), kCounts.data(), kCounts.size(), MPI_INT, MPI_SUM, mComm);
  ///////////////////////////////////////////////////////////////////////////

}
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////  

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int CSRMat::whichSubMatrix(int colID)
{
  for(int i=1; i<mWorldSize; i++)
    {
      // if start col for submatrix i is greater than colID,
      // colID must belong to sumatrix i-1
      if(colID < mSubmatStartCols[i])
      {
        return i-1;
      }
    }

  // if not found in previous submatrices, must be in submatrix mWorldSize-1
  return mWorldSize-1;
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int CSRMat::rowOnWhichProc(int rowID)
{
  for(int i=1; i<mWorldSize; i++)
    {
      // if start col for submatrix i is greater than colID,
      // colID must belong to sumatrix i-1
      if(rowID < mStartRowOnProc[i])
      {
        return i-1;
      }
    }

  // if not found in previous submatrices, must be in submatrix mWorldSize-1
  return mWorldSize-1;
}
//////////////////////////////////////////////////////////////////////////////
