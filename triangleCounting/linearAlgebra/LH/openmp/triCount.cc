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
// File:      miniTri.cc                                                    //
// Project:   miniTri                                                       //
// Author:    Michael Wolf                                                  //
//                                                                          //
// Description:                                                             //
//              Driver for miniTri miniapp                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cstdlib>
#include <cassert>

#include <omp.h>

#include <sys/time.h>

#include "Graph.h"

//////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  struct timeval t1, t2;

 
  if(argc!=4 && argc!=5)
  {
    std::cerr << "Usage: miniTri.exe matrixFile blockSize numThreads [fileformat ={MM || Bin}]" << std::endl;
    exit(1);
  }

  std::string mat = argv[1];
  int blockSize = atoi(argv[2]);
  int numThreads = atoi(argv[3]);
  bool isBinFile = false;

  if(argc==5)
  {
    std::string fileFormat = std::string(argv[4]);

    if(fileFormat == "MM")
    {
      isBinFile=false; 
    }
    else if(fileFormat == "Bin")
    {
      isBinFile=true;
    }
    else
    {
      std::cerr << "File format must be MM or Bin" << std::endl;
      exit(1);
    }
  }

  omp_set_num_threads(numThreads);

  gettimeofday(&t1, NULL);


  Graph g(mat,isBinFile,blockSize);

  g.countTriangles();


  gettimeofday(&t2, NULL);



  std::cout << "Number of Triangles: " << g.getNumTriangles() << std::endl;

  double eTime = t2.tv_sec - t1.tv_sec + ((t2.tv_usec-t1.tv_usec)/1000000.0);

  std::cout << "TIME - Time to count the number of triangles: " << eTime << std::endl;


  //MMW need to unpermute matrix

}
//////////////////////////////////////////////////////////////////////////////

