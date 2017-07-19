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

#include "CSRmatrix.hpp"
#include "mmUtil.h"
#include "binFileReader.h"

int addNZ(std::map<int,std::list<int> > &nzMap,int col, int elemToAdd);

unsigned int choose2(unsigned int k);


//////////////////////////////////////////////////////////////////////////////
// print function -- outputs matrix to file
//                -- accepts optional filename, "CSRmatrix.out" default name
//////////////////////////////////////////////////////////////////////////////
void CSRMat::print() const
{
  std::cout << "Matrix: " << m << " " << n << " " << nnz << std::endl;

  for(int rownum=0; rownum<m; rownum++)
    {
      for(int nzIdx=0; nzIdx<nnzInRow[rownum]; nzIdx++)
	{
	  std::cout << rownum << " " << cols[rownum][nzIdx] << " { ";

	  std::cout << vals[rownum][nzIdx] << std::endl;
	  if(vals2.size()!=0)
	    {
	      std::cout << ", " << vals2[rownum][nzIdx];
	    }
	  std::cout << "}" << std::endl;
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

  for(int rownum=0; rownum<m; rownum++)
    {
      int nrows = nnzInRow[rownum];
      for(int nzIdx=0; nzIdx<nrows; nzIdx++)
	{
	  matList.push_back(rownum);

	  matList.push_back(vals[rownum][nzIdx]);

	  // perhaps should add check for this
	  matList.push_back(vals2[rownum][nzIdx]);

	}
    }

  return matList;
}
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// matmat -- level 3 basic linear algebra subroutine  
//        -- Z = AB where Z = this
//////////////////////////////////////////////////////////////////////////////
void CSRMat::matmat(const CSRMat &A, const CSRMat &B)
{
  //////////////////////////////////////////////////////////
  // set dimensions of matrix, build arrays nnzInRow, vals, cols
  //////////////////////////////////////////////////////////
  int oldM=m;
  m = A.getM();
  n = B.getN();

  if(oldM!=m)
  {
    nnzInRow.resize(m);
    cols.resize(m);
    vals.resize(m);
    vals2.resize(m);
  }
  //////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  // Compute matrix entries one row at a time
  ///////////////////////////////////////////////////////////////////////////
  nnz =0;

  int tmpNNZ =0;

  for (int rownum=0; rownum<m; rownum++)
  {
    nnzInRow[rownum]=0;
    std::map<int,std::list<int> > newNZs;

    int nnzInRowA = A.getNNZInRow(rownum);

    for(int nzindxA=0; nzindxA<nnzInRowA; nzindxA++)
    {
      int colA=A.getCol(rownum, nzindxA);

      int nnzInRowB = B.getNNZInRow(colA);

      for(int nzindxB=0; nzindxB<nnzInRowB; nzindxB++)
      {
        int colB=B.getCol(colA, nzindxB);

        nnzInRow[rownum] += addNZ(newNZs,colB, colA);
      }
    }

    /////////////////////////////////////////////////
    // Strip out any nonzeros that have only one element
    //   This is an optimization for Triangle Enumeration
    //   Algorithm #2                                
    /////////////////////////////////////////////////
    std::map<int,std::list<int> >::iterator iter;

    for (iter=newNZs.begin(); iter!=newNZs.end(); )
    {
      if((*iter).second.size()==1)
      {
        newNZs.erase(iter++); // Remove nonzero
        nnzInRow[rownum]--;   // One less nonzero in row
      }
      else
      {
	++iter;
      }
    }
    /////////////////////////////////////////////////



    // If there are nonzeros in this row
    if(nnzInRow[rownum] > 0)
    {
      /////////////////////////////////////////
      //Allocate memory for this row
      /////////////////////////////////////////
      cols[rownum].resize(nnzInRow[rownum]);
      vals[rownum].resize(nnzInRow[rownum]);
      vals2[rownum].resize(nnzInRow[rownum]);

      /////////////////////////////////////////
      //Copy new data into row
      /////////////////////////////////////////
      std::map<int,std::list<int> >::iterator iter;
      int nzcnt=0;

      // Iterate through list
      for (iter=newNZs.begin(); iter!=newNZs.end(); iter++)
      {
	  cols[rownum][nzcnt]= (*iter).first;

	  std::list<int>::const_iterator lIter=(*iter).second.begin();
	  vals[rownum][nzcnt] = *lIter;
	  lIter++;
	  vals2[rownum][nzcnt]= *lIter;

	  nzcnt++;
      }
      /////////////////////////////////////////
    }
    /////////////////////////////////////////

    tmpNNZ += nnzInRow[rownum]; 


  } // end loop over rows
  ///////////////////////////////////////////////////////////////////////////

  nnz = tmpNNZ;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CSRMat::readMMMatrix(const char *fname)
{
  //////////////////////////////////////////////////////////////
  // Build edge list from MM file                               
  //////////////////////////////////////////////////////////////
  int numVerts;
  int numEdges;
  std::vector<edge_t> edgeList;

  buildEdgeListFromMM(fname, numVerts, numEdges, edgeList);

  m = numVerts;
  n = numVerts;
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Allocate memory for matrix structure                       
  //////////////////////////////////////////////////////////////
  std::vector< std::map<int,int> > rowSets(m);
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Copy data from edgelist to temporary row structure         
  //////////////////////////////////////////////////////////////
  int tmpval=1;

  for (int i=0; i<numEdges; i++)
  {
    if(type==UNDEFINED)
    {
      rowSets[edgeList[i].v0-1][edgeList[i].v1-1] = tmpval;
    }
    else if(type==LOWERTRI)
    {
      if(edgeList[i].v0>edgeList[i].v1)
      {
        rowSets[edgeList[i].v0-1][edgeList[i].v1-1] = tmpval;
      }
    }
    else if(type==UPPERTRI)
    {
      if(edgeList[i].v0<edgeList[i].v1)
      {
        rowSets[edgeList[i].v0-1][edgeList[i].v0-1] = tmpval;
      }
    }
  }
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Free edge lists                                            
  //////////////////////////////////////////////////////////////
  std::vector<edge_t>().swap(edgeList);
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Allocate memory for matrix
  //////////////////////////////////////////////////////////////
  nnzInRow.resize(m);
  cols.resize(m);
  vals.resize(m);
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // copy data from temporary sets to matrix data structures
  //////////////////////////////////////////////////////////////
  std::map<int,int>::iterator iter;

  nnz = 0;
  for(int rownum=0; rownum<m; rownum++)
  {
    int nnzIndx=0;
    int nnzToAdd = rowSets[rownum].size();
    nnzInRow[rownum] = nnzToAdd;
    nnz += nnzToAdd;

    cols[rownum].resize(nnzToAdd);
    vals[rownum].resize(nnzToAdd);

    for (iter=rowSets[rownum].begin();iter!=rowSets[rownum].end();iter++)
    {
	cols[rownum][nnzIndx] = (*iter).first;
	vals[rownum][nnzIndx] = (*iter).second;
	nnzIndx++;
    }

  }
  //////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CSRMat::readBinMatrix(const char *fname)
{
  //////////////////////////////////////////////////////////////
  // Build edge list from MM file                               
  //////////////////////////////////////////////////////////////
  int64_t numVerts=0;
  int64_t numEdges=0;

  std::vector<edge_t> edgeList;

  readBinEdgeFile(fname, numVerts, numEdges, edgeList);

  m = numVerts;
  n = numVerts;
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Allocate memory for matrix structure                       
  //////////////////////////////////////////////////////////////
  std::vector< std::map<int,int> > rowSets(m);
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Copy data from edgelist to temporary row structure         
  //////////////////////////////////////////////////////////////
  int base=0;
  int tmpval=1;

  for (int i=0; i<numEdges; i++)
  {
    if(type==UNDEFINED)
    {
      rowSets[edgeList[i].v0-base][edgeList[i].v1-base] = tmpval;
      rowSets[edgeList[i].v1-base][edgeList[i].v0-base] = tmpval;
    }
    else if(type==LOWERTRI)
    {
      if(edgeList[i].v0>edgeList[i].v1)
      {
        rowSets[edgeList[i].v0-base][edgeList[i].v1-base] = tmpval;
      }
    }
    else if(type==UPPERTRI)
    {
      if(edgeList[i].v0<edgeList[i].v1)
      {
        rowSets[edgeList[i].v0-base][edgeList[i].v0-base] = tmpval;
      }
    }
  }
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Free edge lists                                            
  //////////////////////////////////////////////////////////////
  std::vector<edge_t>().swap(edgeList);
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Allocate memory for matrix
  //////////////////////////////////////////////////////////////
  nnzInRow.resize(m);
  cols.resize(m);
  vals.resize(m);
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // copy data from temporary sets to matrix data structures
  //////////////////////////////////////////////////////////////
  std::map<int,int>::iterator iter;

  nnz = 0;
  for(int rownum=0; rownum<m; rownum++)
  {
    int nnzIndx=0;
    int nnzToAdd = rowSets[rownum].size();
    nnzInRow[rownum] = nnzToAdd;
    nnz += nnzToAdd;

    cols[rownum].resize(nnzToAdd);
    vals[rownum].resize(nnzToAdd);

    for (iter=rowSets[rownum].begin();iter!=rowSets[rownum].end();iter++)
    {
      cols[rownum][nnzIndx] = (*iter).first;
      vals[rownum][nnzIndx] = (*iter).second;
      nnzIndx++;
    }
  }
  //////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CSRMat::createTriMatrix(const CSRMat &matSrc, matrixtype mtype)
{

  m = matSrc.getM();
  n = matSrc.getN();

  assert(m==n);
  assert(mtype==LOWERTRI || mtype==UPPERTRI);
  type = mtype;
  nnz = 0;
  
  //////////////////////////////////////////////////////////////
  // Allocate memory for matrix -- assumes arrays not allocated
  //////////////////////////////////////////////////////////////
  nnzInRow.resize(m);
  cols.resize(m);
  vals.resize(m);
  //////////////////////////////////////////////////////////////

  for(int rownum=0; rownum<m; rownum++)
  {
    std::map<int,int> nzMap;

    int nnzInRowSrc = matSrc.getNNZInRow(rownum);

    for(int nzindxSrc=0; nzindxSrc<nnzInRowSrc; nzindxSrc++)
    {
      int colSrc=matSrc.getCol(rownum, nzindxSrc);
      int valSrc=matSrc.getVal(rownum, nzindxSrc);
         
      // WARNING: assumes there is only 1 element in value for now
      if(type==LOWERTRI && rownum>colSrc)
      {
        nzMap[colSrc]=valSrc;
      }
      else if(type==UPPERTRI && rownum<colSrc)
      {
        nzMap[colSrc]=valSrc;
      }

    }

    int nnzIndx=0;
    int nnzToAdd = nzMap.size();
    nnzInRow[rownum] = nnzToAdd;
    nnz += nnzToAdd;

    cols[rownum].resize(nnzToAdd);
    vals[rownum].resize(nnzToAdd);

    std::map<int,int>::const_iterator iter;

    for (iter=nzMap.begin();iter!=nzMap.end();iter++)
    {
      cols[rownum][nnzIndx] = (*iter).first;
      vals[rownum][nnzIndx] = (*iter).second;
      nnzIndx++;
    }

  } // end of loop over rows

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void CSRMat::createIncidenceMatrix(const CSRMat &matSrc, std::map<int,std::map<int,int> > & eIndices)
{

  m = matSrc.getM();

  assert(type==INCIDENCE);
  nnz = 0;
  
  //////////////////////////////////////////////////////////////
  // Allocate memory for matrix -- assumes arrays not allocated
  //////////////////////////////////////////////////////////////
  nnzInRow.resize(m);
  cols.resize(m);
  vals.resize(m);
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Store columns that need nonzeros
  //////////////////////////////////////////////////////////////
  std::vector<std::set<int> > colsInRow(m);

  int eCnt=0;
  for(int rownum=0; rownum<m; rownum++)
  {

    int nnzInRowSrc = matSrc.getNNZInRow(rownum);

    for(int nzindxSrc=0; nzindxSrc<nnzInRowSrc; nzindxSrc++)
    {
      int colnum=matSrc.getCol(rownum, nzindxSrc);

      if(rownum < colnum)
      {
        colsInRow[rownum].insert(eCnt);
        colsInRow[colnum].insert(eCnt);
	eIndices[rownum][colnum]=eCnt;

        eCnt++;
      }   
      
    }
  }
  n=eCnt;
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // Copy data into matrix data structures
  //     -- Can probably free tmp data structures throughout
  //////////////////////////////////////////////////////////////
  for(int rownum=0; rownum<m; rownum++)
  {
    int nnzToAdd = colsInRow[rownum].size();
    nnzInRow[rownum] = nnzToAdd;
    nnz += nnzToAdd;

    cols[rownum].resize(nnzToAdd);
    vals[rownum].resize(nnzToAdd);

    std::set<int>::const_iterator iter;

    int nnzIndx=0;
    for (iter=colsInRow[rownum].begin();iter!=colsInRow[rownum].end();iter++)
    {
      cols[rownum][nnzIndx] = (*iter);
      vals[rownum][nnzIndx] = 1;
      nnzIndx++;
    }
  }
  //////////////////////////////////////////////////////////////

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

