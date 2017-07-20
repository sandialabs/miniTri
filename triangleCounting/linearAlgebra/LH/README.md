miniTri is a simple, triangle-based data analytics code.  miniTri is a miniapp
in the Mantevo project (http://www.mantevo.org) at Sandia National Laboratories
The primary authors of miniTri are Jon Berry and Michael Wolf (mmwolf@sandia.gov).

miniTri v. 1.0. Copyright (2016) Sandia Corporation.

For questions, contact Jon Berry (jberry@sandia.gov) or Michael Wolf (mmwolf@sandia.gov).

Please read the accompanying README and LICENSE files.

------------------------------------------------
triangleCounting/linearAlgebra/LH:
------------------------------------------------

This directory contains different implementations of a linear algebra based formulation of
triangle counting:

  1. C = (A * H), where A is the adjacency matrix of the graph and H is the incidence matrix of the graph
  2. Number of triangles = the number of entries in C such that C(i,j)=2 divided by 3

A detailed description of this formulation can be found in the following paper:

> Wolf, M.M., J.W. Berry, and D.T. Stark. "A task-based linear algebra Building Blocks 
> approach for scalable graph analytics." High Performance Extreme Computing Conference (HPEC), 
> 2015 IEEE. IEEE, 2015. ([link](http://ieeexplore.ieee.org/document/7322450/?arnumber=7322450))

We have developed several different implementations of this linear algebra-based miniTri using 
different fundamental programming models.  These implementaions are organized as follows:

* __openmp__ -- data parallel OpenMP implementation
* __serial__ -- serial reference implementation 




