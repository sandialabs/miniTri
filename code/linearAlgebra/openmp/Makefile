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
UTILDIRS = ../../utils/
INCDIRS = -I. -I$(UTILDIRS)
CCC = g++ 
BOOSTCCFLAGS = -I/opt/boost-1.55/include
CCFLAGS = -O3 -Wall -DNDEBUG -fopenmp $(BOOSTCCFLAGS)
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

all 		:	lib miniTri
matrix.o	:	CSRMatrix.h 

lib		:	$(LIBOBJECTS) $(UTILOBJECTS)
	ar rvu libSPLA.a $(LIBOBJECTS) $(UTILOBJECTS)

miniTri	:	lib CSRMatrix.h Graph.h
	$(CCC) $(INCDIRS) $(LIBPATH) $(CCFLAGS) -o miniTri.exe miniTri.cc -lSPLA 

clean	:
	rm -f *.o *~ libSPLA.a *.out miniTri.exe
