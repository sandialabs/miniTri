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

int addNZ(std::map<int,std::list<int> > &nzMap,int col, int elemToAdd);

std::list<int> invalidList;

void createPermutation(int *degree, std::vector<int> &perm, std::vector<int> &iperm);
void formDegreeMultiMap(int *degree, int size, std::multimap<int,int> &degreeMap);


//////////////////////////////////////////////////////////////////////////////
// print function -- outputs matrix to file
//                -- accepts optional filename, "CSRmatrix.out" default name
//////////////////////////////////////////////////////////////////////////////
void CSRMat::print() const
{
  std::list<int>::const_iterator it;

  std::cout << "Matrix: " << m << " " << n << " " << nnz << std::endl;

  for(int rownum=0; rownum<m; rownum++)
  {
    for(int nzIdx=0; nzIdx<nnzInRow[rownum]; nzIdx++)
    {
      std::cout << rownum << " " << cols[rownum][nzIdx] << " { ";

      for (it = vals[rownum][nzIdx].begin(); it != vals[rownum][nzIdx].end(); ++it) 
      {
	std::cout << *it << " ";
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
  std::list<int>::const_iterator it;

  for(int rownum=0; rownum<m; rownum++)
  {
    int nrows = nnzInRow[rownum];
    for(int nzIdx=0; nzIdx<nrows; nzIdx++)
    {

      for (it = vals[rownum][nzIdx].begin(); it != vals[rownum][nzIdx].end(); ++it) 
      {
        matList.push_back(*it);
      }
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
    vals = new std::list<int>*[m];
    
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

    // If there are nonzeros in this row
    if(nnzInRow[rownum] > 0)
    {
      /////////////////////////////////////////
      //Allocate memory for this row
      /////////////////////////////////////////

      cols[rownum] = new int[nnzInRow[rownum]];
      vals[rownum] = new std::list<int>[nnzInRow[rownum]];

      /////////////////////////////////////////
      //Copy new data into row
      /////////////////////////////////////////
      std::map<int,std::list<int> >::iterator iter;
      int nzcnt=0;

      // Iterate through list
      for (iter=newNZs.begin(); iter!=newNZs.end(); iter++)
      {
        cols[rownum][nzcnt]= (*iter).first;

        const std::list<int> & lstRef = (*iter).second;

	std::list<int>::const_iterator lIter;
	for (lIter=lstRef.begin(); lIter!=lstRef.end(); lIter++)
	{
          vals[rownum][nzcnt].push_back(*lIter);
	}
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
    std::map<int,std::list<int> > newNZs;
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
	std::list<int> newval;

        const std::list<int> &lstRef = vals[rownum][nzIndx];        
        std::list<int>::const_iterator lstIter;

        // Iterate through list
        // For each element, insert triangles (row,element,col)
        for (lstIter=lstRef.begin(); lstIter!=lstRef.end(); lstIter++)
        {
          newval.push_back(rownum);
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
      vals[rownum] = new std::list<int>[nnzInRow[rownum]];
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
      cols[rownum][nzIndx] = (*mapIt).first;

      // clear value list
      vals[rownum][nzIndx].clear();

      const std::list<int> &lstRef = (*mapIt).second;

      // Iterate through list
      // For each element, insert triangles (row,element,col)
      for (lstIter=lstRef.begin(); lstIter!=lstRef.end(); lstIter++)
      {
        vals[rownum][nzIndx].push_back(*lstIter);
      }

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
  vals = new std::list<int>*[m];
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
    vals[rownum] = new std::list<int> [nnzToAdd];

    for (iter=rowSets[rownum].begin();iter!=rowSets[rownum].end();iter++)
    {
      cols[rownum][nnzIndx] = (*iter).first;
      vals[rownum][nnzIndx].push_back( (*iter).second);
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
  vals = new std::list<int>*[m];
  //////////////////////////////////////////////////////////////

  for(int rownum=0; rownum<m; rownum++)
  {
    std::map<int,int> nzMap;

    int nnzInRowSrc = matSrc.getNNZInRow(rownum);

    for(int nzindxSrc=0; nzindxSrc<nnzInRowSrc; nzindxSrc++)
    {
      int colSrc=matSrc.getCol(rownum, nzindxSrc);
      const std::list<int> & valSrc = matSrc.getVal(rownum, nzindxSrc);
         
      // WARNING: assumes there is only 1 element in value for now
      if(type==LOWERTRI && rownum>colSrc)
      {
        nzMap[colSrc]=valSrc.front();
      }
      else if(type==UPPERTRI && rownum<colSrc)
      {
        nzMap[colSrc]=valSrc.front();
      }

    }

    int nnzIndx=0;
    int nnzToAdd = nzMap.size();
    nnzInRow[rownum] = nnzToAdd;
    nnz += nnzToAdd;

    cols[rownum] = new int[nnzToAdd];
    vals[rownum] = new std::list<int> [nnzToAdd];

    std::map<int,int>::const_iterator iter;

    for (iter=nzMap.begin();iter!=nzMap.end();iter++)
    {
      cols[rownum][nnzIndx] = (*iter).first;
      vals[rownum][nnzIndx].push_back( (*iter).second);
      nnzIndx++;
    }

  } // end of loop over rows

}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//
// perhaps could improve algorithm by not requiring copy of data
////////////////////////////////////////////////////////////////////////////////
void CSRMat::permute()
{
  std::vector<int> perm(m);
  std::vector<int> iperm(m);

  createPermutation(nnzInRow, perm, iperm);

  ///////////////////////////////////////////////////////////////////////////
  // Set temp pointers to save original order
  ///////////////////////////////////////////////////////////////////////////
  std::vector<int> tmpNNZ(m);
  std::vector<int *> tmpCols(m);
  std::vector<std::list<int> *> tmpVals(m);

  for(int rownum=0; rownum<m; rownum++)
  {
    tmpNNZ[rownum] = nnzInRow[rownum];
    tmpCols[rownum] = cols[rownum];
    tmpVals[rownum] = vals[rownum];
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  //Permute matrix, row by row
  ///////////////////////////////////////////////////////////////////////////
  for(int rownum=0; rownum<m; rownum++)
  {
    //////////////////////////////////////////////////////////////////////
    // Copy data into permuted order
    //////////////////////////////////////////////////////////////////////
    nnzInRow[rownum] = tmpNNZ[iperm[rownum]];    
    cols[rownum] = tmpCols[iperm[rownum]];
    
    /////////////////////////////////////////////////////////////////
    // permute column numbers as well -- perhaps should sort this
    /////////////////////////////////////////////////////////////////
    for(int i=0;i<nnzInRow[rownum];i++)
    {
      cols[rownum][i] = perm[cols[rownum][i]];
    }
    /////////////////////////////////////////////////////////////////

    vals[rownum] = tmpVals[iperm[rownum]];

    //////////////////////////////////////////////////////////////////////

  }
  ///////////////////////////////////////////////////////////////////////////

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


