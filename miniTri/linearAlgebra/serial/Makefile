##############################################################################
#                                                                            #
# File:      Makefile                                                        #
# Project:   miniTri                                                         #
# Author:    Michael Wolf                                                    #
#                                                                            #
# Description:                                                               #
#              Makefile for code.                                            #
#                                                                            #
##############################################################################

UTILDIR = ../../utils/
INCDIRS = -I. -I$(UTILDIR)
CCC = g++ 
CCFLAGS = -O3 -Wall -DNDEBUG -std=c++11
LIBPATH = -L. 

#--------------------------------------------------
LIBSOURCES =             \
          CSRmatrix.cpp   \
          Graph.cpp   

LIBOBJECTS         = $(LIBSOURCES:.cpp=.o) 
UTILOBJECTS	   = mmio.o mmUtil.o binFileReader.o

#--------------------------------------------------

.SUFFIXES:

%.o : %.cc
	$(CCC) $< -c $(CCFLAGS) $(INCDIRS)

%.o : %.cpp
	$(CCC) $< -c $(CCFLAGS) $(INCDIRS)



#--------------------------------------------------
vpath %.cc $(UTILDIR)
vpath %.cpp $(UTILDIR)

all 		:	lib miniTri 
matrix.o	:	CSRmatrix.h 

lib		:	 $(LIBOBJECTS) $(UTILOBJECTS)
	ar rvu libSPLA.a $(LIBOBJECTS) $(UTILOBJECTS)

miniTri	:	lib CSRmatrix.hpp Graph.hpp
	$(CCC) $(INCDIRS) $(LIBPATH) $(CCFLAGS) -o miniTri.exe miniTri.cpp -lSPLA 

clean	:
	rm -f *.o *~ libSPLA.a *.out miniTri.exe
