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

int addNZ(std::map<int,int> &nzMap,int col, int elemToAdd);

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
      std::cout << rownum << " " << cols[rownum][nzIdx] << vals[rownum][nzIdx] << std::endl;;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Sums matrix elements
////////////////////////////////////////////////////////////////////////////////
int CSRMat::getSumElements() const
{
  int sum = 0;

  for(int rownum=0; rownum<m; rownum++)
  {
    int nrows = nnzInRow[rownum];
    for(int nzIdx=0; nzIdx<nrows; nzIdx++)
    {
      sum +=  vals[rownum][nzIdx];
    }
  }

  return sum;
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
    if(nnzInRow!=0)
    {
      delete [] nnzInRow;
    }

    if(vals!=0)
    {
      for(int i=0; i<oldM; i++)
      {
        if(vals[i]!=0)
        {
          delete [] vals[i];
        }
      }
      delete [] vals;
    }

    if(cols!=0)
    {
      for(int i=0; i<oldM; i++)
      {
        if(cols[i]!=0)
        {
          delete [] cols[i];
        }
      }
      delete [] cols;
    }
  
    nnzInRow = new int[m];
    cols = new int*[m];
    vals = new int*[m];
    
    for(int i=0;i<m;i++)
    {
      cols[i]=0;
      vals[i]=0;
    }


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
    std::map<int,int> newNZs;

    int nnzInRowA = A.getNNZInRow(rownum);

    for(int nzindxA=0; nzindxA<nnzInRowA; nzindxA++)
    {
      int colA=A.getCol(rownum, nzindxA);

      int nnzInRowB = B.getNNZInRow(colA);

      for(int nzindxB=0; nzindxB<nnzInRowB; nzindxB++)
      {
        int colB=B.getCol(colA, nzindxB);

        // For now set val to 1, should really do multiplication
        nnzInRow[rownum] += addNZ(newNZs,colB, 1);
      }
    }

    // If there are nonzeros in this row
    if(nnzInRow[rownum] > 0)
    {
      /////////////////////////////////////////
      //Allocate memory for this row
      /////////////////////////////////////////

      cols[rownum] = new int[nnzInRow[rownum]];
      vals[rownum] = new int[nnzInRow[rownum]];

      /////////////////////////////////////////
      //Copy new data into row
      /////////////////////////////////////////
      std::map<int,int>::iterator iter;
      int nzcnt=0;

      for (iter=newNZs.begin(); iter!=newNZs.end(); iter++)
      {
        cols[rownum][nzcnt]= (*iter).first;
        vals[rownum][nzcnt]= (*iter).second;

        nzcnt++;
      }
      /////////////////////////////////////////
    }
    else
    {
      cols[rownum] = 0;
      vals[rownum] = 0;
    }
    /////////////////////////////////////////

    tmpNNZ += nnzInRow[rownum]; 


  } // end loop over rows
  ///////////////////////////////////////////////////////////////////////////

  nnz = tmpNNZ;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// EWMult -- A = A .* W, where A = this
////////////////////////////////////////////////////////////////////////////////
void CSRMat::EWMult(const CSRMat &W)
{
  //////////////////////////////////////////////////////////
  // check dimensions of matrices
  //////////////////////////////////////////////////////////
  assert(m == W.getM());
  assert(n == W.getN());
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  // Calculate matrix values for future matrix
  //////////////////////////////////////////////////////////
  nnz =0;
  // Loop over each row
  for(int rownum=0;rownum<m;rownum++)
  {
    std::map<int,int> newNZs;
    std::set<int> colsWrow;

    ////////////////////////////////////////////////
    // Insert cols from W row into temp data structure colsWrow
    ////////////////////////////////////////////////    
    int nnzInRowW = W.getNNZInRow(rownum);

    for(int nzIndxW=0; nzIndxW<nnzInRowW; nzIndxW++)
    {
      int colW = W.getCol(rownum, nzIndxW);

      colsWrow.insert(colW);
    }
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // Loop of nzs in row for this matrix.
    //     If column matches column in W (in colsWrow),
    //     this results in nonzero in new matrix
    ////////////////////////////////////////////////
    std::set<int>::const_iterator it;
    for (int nzIndx=0; nzIndx<nnzInRow[rownum]; nzIndx++)
    {
      int colIndx = cols[rownum][nzIndx];

      it=colsWrow.find(colIndx);      

      //////////////////////////////////////
      //If columns match, multiply nonzero elements
      //////////////////////////////////////      
      if(it != colsWrow.end())
      {
        newNZs.insert(std::pair<int,int>(colIndx,vals[rownum][nzIndx]));
      }

    }
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // Rebuild current matrix data structures with new data
    ////////////////////////////////////////////////
    int oldNNZ = nnzInRow[rownum];
    nnzInRow[rownum] = newNZs.size();  
    nnz += nnzInRow[rownum];
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // Rellocate memory for col and val arrays if necessary
    ////////////////////////////////////////////////
    if(oldNNZ != nnzInRow[rownum])
    {
      if(cols[rownum]!=0)
      {
        delete [] cols[rownum];
      }
      cols[rownum] = new int[nnzInRow[rownum]];

      if(vals[rownum]!=0)
      {
        delete [] vals[rownum];
      }
      vals[rownum] = new int[nnzInRow[rownum]];
    }
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // Copy new nz data into matrix
    ////////////////////////////////////////////////
    int nzIndx=0;

    std::map<int,int>::const_iterator mapIt;
    std::list<int>::const_iterator lstIter;

    for(mapIt=newNZs.begin();mapIt!=newNZs.end();++mapIt)
    {
      cols[rownum][nzIndx] = (*mapIt).first;
      vals[rownum][nzIndx] = (*mapIt).second;

      nzIndx++;
    }
    //////////////////////////////////////////////////
  }
  //////////////////////////////////////////////////////////

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
  nnz = numEdges;
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

  for (int i=0; i<nnz; i++)
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
  nnzInRow = new int[m];
  cols = new int*[m];
  vals = new int*[m];
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

    cols[rownum] = new int[nnzToAdd];
    vals[rownum] = new int[nnzToAdd];

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
  nnzInRow = new int[m];
  cols = new int*[m];
  vals = new int*[m];
  //////////////////////////////////////////////////////////////

  for(int rownum=0; rownum<m; rownum++)
  {
    std::map<int,int> nzMap;

    int nnzInRowSrc = matSrc.getNNZInRow(rownum);

    for(int nzindxSrc=0; nzindxSrc<nnzInRowSrc; nzindxSrc++)
    {
      int colSrc = matSrc.getCol(rownum, nzindxSrc);
      int valSrc = matSrc.getVal(rownum, nzindxSrc);
         
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

    cols[rownum] = new int[nnzToAdd];
    vals[rownum] = new int[nnzToAdd];

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


//////////////////////////////////////////////////////////////////////////////
// addNZ -- For a given row, add a column for a nonzero into a sorted list
//////////////////////////////////////////////////////////////////////////////
int addNZ(std::map<int,int> &nzMap,int col, int elemToAdd)
{
  std::map<int,int>::iterator it;

  it = nzMap.find(col);

  //////////////////////////////////////
  //If columns match, no additional nz, add element to end of list
  //////////////////////////////////////      
  if(it != nzMap.end())
  {
    (*it).second += elemToAdd;
    return 0;
  }

  nzMap.insert(std::pair<int,int>(col, elemToAdd));
  return 1;
}
//////////////////////////////////////////////////////////////////////////////


