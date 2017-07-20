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
UTILDIRS = ../../../utils/
INCDIRS = -I. -I$(UTILDIRS)
CCC = g++ 
CCFLAGS = -O3 -Wall -DNDEBUG -fopenmp -std=c++11
LIBPATH = -L. 


#--------------------------------------------------
LIBSOURCES =             \
          CSRMatrix.cc   \
          Graph.cc 

LIBOBJECTS         = $(LIBSOURCES:.cc=.o) 
UTILOBJECTS        = mmio.o mmUtil.o binFileReader.o

#--------------------------------------------------

.SUFFIXES:

%.o : %.cc
	$(CCC) $< -c $(CCFLAGS) $(INCDIRS)

%.o : %.cpp
	$(CCC) $< -c $(CCFLAGS) $(INCDIRS)

%.o : %.f
	$(F77) $< -c $(CCFLAGS) $(INCDIRS)


#--------------------------------------------------
vpath %.cc $(UTILDIRS)
vpath %.cpp $(UTILDIRS)

all 		:	lib triCount
matrix.o	:	CSRMatrix.h 

lib		:	$(LIBOBJECTS) $(UTILOBJECTS)
	ar rvu libSPLA.a $(LIBOBJECTS) $(UTILOBJECTS)

triCount	:	lib CSRMatrix.h Graph.h
	$(CCC) $(INCDIRS) $(LIBPATH) $(CCFLAGS) -o triCount.exe triCount.cc -lSPLA 

clean	:
	rm -f *.o *~ libSPLA.a *.out miniTri.exe
