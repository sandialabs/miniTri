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

#include "CSRmatrix.h"
#include "mmUtil.h"
#include "mmio.h"

void printSubmat(const CSRSubmat &submat, int startRow,int locNumRows);

int addNZ(std::map<int,std::list<int> > &nzMap,int col, int elemToAdd);

void createPermutation(int *degree, std::vector<int> &perm, std::vector<int> &iperm);
void formDegreeMultiMap(int *degree, int size, std::multimap<int,int> &degreeMap);


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

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat &submat = mSubmat[submatNum];

    printSubmat(submat,mStartRow,mLocNumRows);

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

      for(int nzIdx=0; nzIdx<submat.nnzInRow[rownum]; nzIdx++)
      {

        std::cout << startRow + rownum << " " << submat.cols[rownum][nzIdx] << " { ";


        for (it = submat.vals[rownum][nzIdx].begin(); it != submat.vals[rownum][nzIdx].end(); ++it) 
        {
	  std::cout << *it << " ";
        }
        std::cout << "}" << std::endl;
      }
  }

}
//////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Sums matrix elements -- Each processor calculates their share, no reduction.
////////////////////////////////////////////////////////////////////////////////
std::list<int> CSRMat::getSumElements() const
{
  std::list<int> matList;
  std::list<int>::const_iterator it;

  ///////////////////////////////////////////////////////////////////////////
  // Local operation
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat &submat = mSubmat[submatNum];

    for(int rownum=0; rownum<mLocNumRows; rownum++)
    {
      int nrows = submat.nnzInRow[rownum];
      for(int nzIdx=0; nzIdx<nrows; nzIdx++)
      {
        for (it = submat.vals[rownum][nzIdx].begin(); it != submat.vals[rownum][nzIdx].end(); ++it) 
        {
          matList.push_back(*it);
        }
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  return matList;
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

  if(oldGlobNumRows!=mGlobNumRows)
  {
    for(unsigned int submatNum=0; submatNum< mSubmat.size(); submatNum++)
    {
      freeSubmat(mSubmat[submatNum],oldGlobNumRows);      
    }

    mSubmat.resize(mWorldSize);

    for(unsigned int submatNum=0; submatNum< mSubmat.size(); submatNum++)
    {
      mSubmat[submatNum].nnzInRow = new int[mLocNumRows];
      mSubmat[submatNum].cols = new int*[mLocNumRows];
      mSubmat[submatNum].vals = new std::list<int>*[mLocNumRows];
    
      for(int i=0;i<mLocNumRows;i++)
      {
        mSubmat[submatNum].cols[i]=0;
        mSubmat[submatNum].vals[i]=0;
      }

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

      freeSubmat(remSubmatB,remNumRowsB);
      //////////////////////////////////////////////////////////////////////
    }

  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Copy resulting nonzeros into matrix object
  ///////////////////////////////////////////////////////////////////////////
  mLocNNZ=0;
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    const CSRSubmat & submatC = mSubmat[submatNum];

    for (int rownum=0; rownum<mLocNumRows; rownum++)
    {
      const std::map<int,std::list<int> > &nzMap = CNZs[submatNum][rownum];


      submatC.nnzInRow[rownum] = nzMap.size();

      // If there are nonzeros in this row
      if(submatC.nnzInRow[rownum] > 0)
      {
        /////////////////////////////////////////
        //Allocate memory for this row
        /////////////////////////////////////////
        submatC.cols[rownum] = new int[submatC.nnzInRow[rownum]];
        submatC.vals[rownum] = new std::list<int>[submatC.nnzInRow[rownum]];

        /////////////////////////////////////////
        //Copy new data into row
        /////////////////////////////////////////
        std::map<int,std::list<int> >::const_iterator iter;
        int nzcnt=0;

        // Iterate through map
        for (iter=nzMap.begin(); iter!=nzMap.end(); iter++)
        {
          submatC.cols[rownum][nzcnt]= (*iter).first;

          const std::list<int> & lstRef = (*iter).second;

 	  std::list<int>::const_iterator lIter;
	  for (lIter=lstRef.begin(); lIter!=lstRef.end(); lIter++)
	  {
            submatC.vals[rownum][nzcnt].push_back(*lIter);
	  }
          nzcnt++;
        }
        /////////////////////////////////////////
      }
      else
      {
        submatC.cols[rownum] = 0;
        submatC.vals[rownum] = 0;
      }
      /////////////////////////////////////////

      mLocNNZ += submatC.nnzInRow[rownum];  

    } // Loop over rows

  } // loop over submats
  ///////////////////////////////////////////////////////////////////////////

  MPI_Allreduce(&mLocNNZ, &mGlobNNZ, 1, MPI_INT, MPI_SUM,mComm);

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// EWMult -- A = A .* W, where A = this
//        -- Completely local operation, except for MPI_Allreduce to sum up NNZ
////////////////////////////////////////////////////////////////////////////////
void CSRMat::EWMult(const CSRMat &W)
{
  //////////////////////////////////////////////////////////////
  // check dimensions of matrices
  //////////////////////////////////////////////////////////////
  assert(mGlobNumRows == W.getGlobNumRows());
  assert(mGlobNumCols == W.getGlobNumCols());
  assert(mLocNumRows == W.getLocNumRows());
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Calculate matrix values for future matrix
  //////////////////////////////////////////////////////////////
  mLocNNZ =0;

  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    CSRSubmat &submat = mSubmat[submatNum];

    // Loop over each row
    for(int rownum=0;rownum<mLocNumRows;rownum++)
    {
      std::map<int,std::list<int> > newNZs;
      std::set<int> colsWrow;

      ////////////////////////////////////////////////
      // Insert cols from W row into temp data structure colsWrow
      ////////////////////////////////////////////////    
      int nnzInRowW = W.getNNZInRow(submatNum,rownum);

      for(int nzIndxW=0; nzIndxW<nnzInRowW; nzIndxW++)
      {
        int colW = W.getCol(submatNum,rownum, nzIndxW);

        colsWrow.insert(colW);
      }
      ////////////////////////////////////////////////

      ////////////////////////////////////////////////
      // Loop of nzs in row for this matrix.
      //     If column matches column in W (in colsWrow),
      //     this results in nonzero in new matrix
      ////////////////////////////////////////////////
      std::set<int>::const_iterator it;
      for (int nzIndx=0; nzIndx<submat.nnzInRow[rownum]; nzIndx++)
      {
        int colIndx = submat.cols[rownum][nzIndx];

        it=colsWrow.find(colIndx);      

        //////////////////////////////////////
        //If columns match, multiply nonzero elements
        //////////////////////////////////////      
        if(it != colsWrow.end())
        {
	  std::list<int> newval;

          const std::list<int> &lstRef = submat.vals[rownum][nzIndx];        
          std::list<int>::const_iterator lstIter;

          // Iterate through list
          // For each element, insert triangles (row,element,col)
          for (lstIter=lstRef.begin(); lstIter!=lstRef.end(); lstIter++)
          {
            newval.push_back(mStartRow+rownum);
            newval.push_back(*lstIter);
            newval.push_back(colIndx);
          }

          //MMW: does this copy properly
          newNZs.insert(std::pair<int,std::list<int> >(colIndx,newval));
        }
      }
      ////////////////////////////////////////////////

      ////////////////////////////////////////////////
      // Rebuild current matrix data structures with new data
      ////////////////////////////////////////////////
      int oldNNZ = submat.nnzInRow[rownum];
      submat.nnzInRow[rownum] = newNZs.size();  
      mLocNNZ += submat.nnzInRow[rownum];
      ////////////////////////////////////////////////

      ////////////////////////////////////////////////
      // Rellocate memory for col and val arrays if necessary
      ////////////////////////////////////////////////
      if(oldNNZ != submat.nnzInRow[rownum])
      {
        if(submat.cols[rownum]!=0)
        {
          delete [] submat.cols[rownum];
        }
        submat.cols[rownum] = new int[submat.nnzInRow[rownum]];

        if(submat.vals[rownum]!=0)
        {
          delete [] submat.vals[rownum];
        }
        submat.vals[rownum] = new std::list<int>[submat.nnzInRow[rownum]];
      }
      ////////////////////////////////////////////////

      ////////////////////////////////////////////////
      // Copy new nz data into matrix
      ////////////////////////////////////////////////
      int nzIndx=0;

      std::map<int,std::list<int> >::const_iterator mapIt;
      std::list<int>::const_iterator lstIter;

      for(mapIt=newNZs.begin();mapIt!=newNZs.end();++mapIt)
      {
        submat.cols[rownum][nzIndx] = (*mapIt).first;

        // clear value list
        submat.vals[rownum][nzIndx].clear();

        const std::list<int> &lstRef = (*mapIt).second;

        // Iterate through list
        // For each element, insert triangles (row,element,col)
        for (lstIter=lstRef.begin(); lstIter!=lstRef.end(); lstIter++)
        {
          submat.vals[rownum][nzIndx].push_back(*lstIter);
        }

        nzIndx++;

      }
      //////////////////////////////////////////////////
    } //for rownum
    //////////////////////////////////////////////////////////////
  } // for submatnum
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
void CSRMat::readMMMatrix(const char *fname)
{

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
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
  // Sets values for mSubmatStartCols
  ///////////////////////////////////////////////////////////////////////////
  for(int rank=0; rank<mWorldSize; rank++)
  {
    int numcols;
    partitionMatrix(mGlobNumCols,mWorldSize,rank,numcols,mSubmatStartCols[rank]);
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Allocate memory for matrix
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    mSubmat[submatNum].nnzInRow = new int[mLocNumRows];
    mSubmat[submatNum].cols = new int*[mLocNumRows];
    mSubmat[submatNum].vals = new std::list<int>*[mLocNumRows];
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

    if(type==UNDEFINED)
    {
      int subMat = whichSubMatrix(tmpcol-1);
      rowSets[subMat][tmprow-1-mStartRow][tmpcol-1] = tmpval;
    }
    else if(type==LOWERTRI)
    {
      if(tmprow>tmpcol)
      {
        int subMat = whichSubMatrix(tmpcol-1);
        rowSets[subMat][tmprow-1-mStartRow][tmpcol-1] = tmpval;       
      }
    }
    else if(type==UPPERTRI)
    {
      if(tmprow<tmpcol)
      {
        int subMat = whichSubMatrix(tmpcol-1);
        rowSets[subMat][tmprow-1-mStartRow][tmpcol-1] = tmpval;       
      }
    }

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
      submat.nnzInRow[rownum] = nnzToAdd;
      mLocNNZ += nnzToAdd;

      if(nnzToAdd>0)
      {
        submat.cols[rownum] = new int[nnzToAdd];
        submat.vals[rownum] = new std::list<int> [nnzToAdd];

        for (iter=rowSetsRef[rownum].begin();iter!=rowSetsRef[rownum].end();iter++)
        {
          submat.cols[rownum][nnzIndx] = (*iter).first;
          submat.vals[rownum][nnzIndx].push_back( (*iter).second);
          nnzIndx++;
        }
      }
      else
      {
	submat.cols[rownum] = 0;
        submat.vals[rownum] = 0;
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
void CSRMat::createTriMatrix(const CSRMat &matSrc, matrixtype mtype)
{

  mGlobNumRows = matSrc.getGlobNumRows();
  mGlobNumCols = matSrc.getGlobNumCols();
  mLocNumRows = matSrc.getLocNumRows();
  mStartRow = matSrc.getStartRow();

  assert(mGlobNumRows==mGlobNumCols);
  assert(mtype==LOWERTRI || mtype==UPPERTRI);
  type = mtype;
  mLocNNZ = 0;

  ///////////////////////////////////////////////////////////////////////////
  // Allocate memory for matrix -- assumes arrays not allocated
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    mSubmat[submatNum].nnzInRow = new int[mLocNumRows];
    mSubmat[submatNum].cols = new int*[mLocNumRows];
    mSubmat[submatNum].vals = new std::list<int>*[mLocNumRows];
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  for(int submatNum=0; submatNum<mWorldSize; submatNum++)
  {
    for(int rownum=0; rownum<mLocNumRows; rownum++)
    {
      std::map<int,int> nzMap;

      int nnzInRowSrc = matSrc.getNNZInRow(submatNum,rownum);

      for(int nzindxSrc=0; nzindxSrc<nnzInRowSrc; nzindxSrc++)
      {
        int colSrc=matSrc.getCol(submatNum, rownum, nzindxSrc);
        const std::list<int> & valSrc = matSrc.getVal(submatNum, rownum, nzindxSrc);
         
        // WARNING: assumes there is only 1 element in value for now
        if(type==LOWERTRI && rownum+mStartRow>colSrc)
        {
          nzMap[colSrc]=valSrc.front();
        }
        else if(type==UPPERTRI && rownum+mStartRow<colSrc)
        {
          nzMap[colSrc]=valSrc.front();
        }
      }

      int nnzIndx=0;
      int nnzToAdd = nzMap.size();
      mSubmat[submatNum].nnzInRow[rownum] = nnzToAdd;
      mLocNNZ += nnzToAdd;

      if(nnzToAdd>0)
      {
        mSubmat[submatNum].cols[rownum] = new int[nnzToAdd];
        mSubmat[submatNum].vals[rownum] = new std::list<int> [nnzToAdd];

        std::map<int,int>::const_iterator iter;

        for (iter=nzMap.begin();iter!=nzMap.end();iter++)
        {
          mSubmat[submatNum].cols[rownum][nnzIndx] = (*iter).first;
          mSubmat[submatNum].vals[rownum][nnzIndx].push_back( (*iter).second);
          nnzIndx++;
        }
      }
      else
      {
        mSubmat[submatNum].cols[rownum] = 0;
        mSubmat[submatNum].vals[rownum] = 0;
      }

    } // end of loop over rows

  } // end loop over submatrices
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Sum mLocNNZ over all processes to obtain mGlobNNZ
  ///////////////////////////////////////////////////////////////////////////
  MPI_Allreduce(&mLocNNZ,&mGlobNNZ,1,MPI_INT, MPI_SUM, mComm);
  ///////////////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//
// perhaps could improve algorithm by not requiring copy of data
//
// MMW: rework for MPI
////////////////////////////////////////////////////////////////////////////////
void CSRMat::permute()
{
//   std::vector<int> perm(m);
//   std::vector<int> iperm(m);

//   createPermutation(nnzInRow, perm, iperm);

//   ///////////////////////////////////////////////////////////////////////////
//   // Set temp pointers to save original order
//   ///////////////////////////////////////////////////////////////////////////
//   std::vector<int> tmpNNZ(m);
//   std::vector<int *> tmpCols(m);
//   std::vector<std::list<int> *> tmpVals(m);

//   for(int rownum=0; rownum<m; rownum++)
//   {
//     tmpNNZ[rownum] = nnzInRow[rownum];
//     tmpCols[rownum] = cols[rownum];
//     tmpVals[rownum] = vals[rownum];
//   }
//   ///////////////////////////////////////////////////////////////////////////

//   ///////////////////////////////////////////////////////////////////////////
//   //Permute matrix, row by row -- this might hurt data locality (ptr swapping)
//   ///////////////////////////////////////////////////////////////////////////
//   for(int rownum=0; rownum<m; rownum++)
//   {
//     //////////////////////////////////////////////////////////////////////
//     // Copy data into permuted order
//     //////////////////////////////////////////////////////////////////////
//     nnzInRow[rownum] = tmpNNZ[iperm[rownum]];    
//     cols[rownum] = tmpCols[iperm[rownum]];
    
//     /////////////////////////////////////////////////////////////////
//     // permute column numbers as well -- perhaps should sort this
//     /////////////////////////////////////////////////////////////////
//     for(int i=0;i<nnzInRow[rownum];i++)
//     {
//       cols[rownum][i] = perm[cols[rownum][i]];
//     }
//     /////////////////////////////////////////////////////////////////

//     vals[rownum] = tmpVals[iperm[rownum]];

//     //////////////////////////////////////////////////////////////////////

//   }
//   ///////////////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void createPermutation(int *degree, std::vector<int> &perm, std::vector<int> &iperm)
{
   std::multimap<int,int> degreeMMap;
   formDegreeMultiMap(degree,perm.size(),degreeMMap);

   std::multimap<int,int>::const_iterator iter;
   int cnt=0;
   for(iter=degreeMMap.begin(); iter!=degreeMMap.end(); ++iter)
   {
       perm[(*iter).second] = cnt;
       iperm[cnt] = (*iter).second;
       cnt++;
   }
}
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Form sorted multimap of (degree,rownum) pairs, currently sorted in increasing order
//////////////////////////////////////////////////////////////////////////////
void formDegreeMultiMap(int *degree, int size, std::multimap<int,int> &degreeMMap)
{
  for(int i=0; i<size; i++)
  {
    degreeMMap.insert(std::pair<int, int>(degree[i], i));
  }
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Multiplies two submatrices together
// Inserts resulting nonzeros into vector of <col,val> value maps
//////////////////////////////////////////////////////////////////////////////
void serialSubmatrixMult(const CSRSubmat & submatA, int numARows, 
			 const CSRSubmat & submatB, int startRowB, 
			 std::vector<std::map<int,std::list<int> > > & submatCNZs)
{

  ///////////////////////////////////////////////////////////////////////////
  // Compute matrix entries one row at a time
  ///////////////////////////////////////////////////////////////////////////
  for (int rownum=0; rownum<numARows; rownum++)
  {
    int nnzInRowA = submatA.nnzInRow[rownum];

    for(int nzindxA=0; nzindxA<nnzInRowA; nzindxA++)
    {
      int colA = submatA.cols[rownum][nzindxA];

      int nnzInRowB = submatB.nnzInRow[colA-startRowB];

      for(int nzindxB=0; nzindxB<nnzInRowB; nzindxB++)
      {
        int colB=submatB.cols[colA-startRowB][nzindxB];

        addNZ(submatCNZs[rownum], colB, colA);
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Send/recv submatrix
//////////////////////////////////////////////////////////////////////////////
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
  remSubmat.nnzInRow = new int[numRowsRecv];
  remSubmat.cols = new int*[numRowsRecv];
  remSubmat.vals = new std::list<int>*[numRowsRecv];
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Sendrecv nnzInRow data
  /////////////////////////////////////////////////////////////////////////
  MPI_Sendrecv(submatToSend.nnzInRow, numRowsSend, MPI_INT, dst, 0,
               remSubmat.nnzInRow, numRowsRecv, MPI_INT, src, 0,
	       comm, &status);
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // loop over nnzInRow data to allocate cols, vals
  /////////////////////////////////////////////////////////////////////////
  for(int rownum=0; rownum<numRowsRecv; rownum++)
  {
    if(remSubmat.nnzInRow[rownum] > 0)
    {
      remSubmat.cols[rownum] = new int[remSubmat.nnzInRow[rownum]];
      remSubmat.vals[rownum] = new std::list<int>[remSubmat.nnzInRow[rownum]];
    }
    else
    {
      remSubmat.cols[rownum] = 0;
      remSubmat.vals[rownum] = 0;
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
    if(remSubmat.nnzInRow[rownum] > 0)
    {
      MPI_Irecv((remSubmat.cols[rownum]),remSubmat.nnzInRow[rownum], MPI_INT, 
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
    if(submatToSend.nnzInRow[rownum] > 0)
    {
      MPI_Send(submatToSend.cols[rownum], submatToSend.nnzInRow[rownum], 
	       MPI_INT, dst, 0, comm);
    }
  }

  MPI_Waitall(numRowsRecv,requests,stats);
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Fill linear array of values (each value can consist of multiple elements)
  // sendValPtr will point to each beginning of each value
  /////////////////////////////////////////////////////////////////////////

  std::vector<int> sendValPtr;
  int curIdx=0;
  sendValPtr.push_back(curIdx);
  std::vector<int> sendValsPacked;

  for(int rownum=0; rownum<numRowsSend; rownum++)
  {
    int nnzInRow = submatToSend.nnzInRow[rownum];

    for(int idx=0; idx<nnzInRow; idx++)
    {
      const std::list<int> &valRef = submatToSend.vals[rownum][idx];
      curIdx += valRef.size();
      sendValPtr.push_back(curIdx);

      std::list<int>::const_iterator iter;
      for(iter=valRef.begin();iter!=valRef.end();iter++)
      {
        sendValsPacked.push_back(*iter);
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Sendrecv valPtr array
  /////////////////////////////////////////////////////////////////////////
  int sendSize = sendValPtr.size();
  int recvSize;

  MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0,
               &recvSize, 1, MPI_INT, src, 0, comm, &status);

  std::vector<int> recvValPtr(recvSize);

  MPI_Sendrecv(&(sendValPtr[0]), sendSize, MPI_INT, dst, 0,
               &(recvValPtr[0]), recvSize, MPI_INT, src, 0, comm, &status);

  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Sendrecv valsPacked array
  /////////////////////////////////////////////////////////////////////////
  sendSize = sendValsPacked.size();

  MPI_Sendrecv(&sendSize, 1, MPI_INT, dst, 0,
               &recvSize, 1, MPI_INT, src, 0, comm, &status);

  std::vector<int> recvValsPacked(recvSize);

  if(sendSize>0)
  {
    if(recvSize>0)
    {
      MPI_Sendrecv(&(sendValsPacked[0]), sendSize, MPI_INT, dst, 0,
                   &(recvValsPacked[0]), recvSize, MPI_INT, src, 0, comm, &status);
    }
    else
    {
      MPI_Send(&(sendValsPacked[0]), sendSize, MPI_INT, dst, 0, comm);
    }
  }
  else if(recvSize>0)
  {
    MPI_Recv(&(recvValsPacked[0]), recvSize, MPI_INT, src, 0, comm, &status);
  }
  /////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////
  // Unpack values and store in submat structure
  /////////////////////////////////////////////////////////////////////////
  if(recvSize>0)
  {
    int nzIdx=0;

    for(int rownum=0; rownum<numRowsRecv; rownum++)
    {
      int nnzInRow = remSubmat.nnzInRow[rownum];

      for(int i=0; i<nnzInRow; i++)
      {
        for(int valElemIdx=recvValPtr[nzIdx]; valElemIdx<recvValPtr[nzIdx+1]; valElemIdx++)
        {
          remSubmat.vals[rownum][i].push_back(recvValsPacked[valElemIdx]);  
        }

        nzIdx++;
      }
    }
  }
  /////////////////////////////////////////////////////////////////////////

}
//////////////////////////////////////////////////////////////////////////////


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

//////////////////////////////////////////////////////////////////////////////
// Frees submatrix
//////////////////////////////////////////////////////////////////////////////
void freeSubmat(CSRSubmat &submat, int nrows)
{
    if(submat.nnzInRow!=0)
    {
      delete [] submat.nnzInRow; submat.nnzInRow=0;
    }

    if(submat.cols!=0)
    {
      for(int i=0;i<nrows;i++)
      {
        if(submat.cols[i]!=0)
	{
          delete [] submat.cols[i]; submat.cols[i]=0;
	}
      }
      delete [] submat.cols; submat.cols=0;
    }

    if(submat.vals!=0)
    {
      for(int i=0;i<nrows;i++)
      {
        if(submat.vals[i]!=0)
	{
          delete [] submat.vals[i]; submat.vals[i]=0;
	}
      }
      delete [] submat.vals; submat.vals=0;
    }

}
//////////////////////////////////////////////////////////////////////////////

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
