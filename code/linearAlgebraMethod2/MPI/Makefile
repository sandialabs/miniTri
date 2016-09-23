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
UTILDIR = ../../../utils/
INCDIRS = -I. -I$(UTILDIR)
CCC = mpicxx 
CCFLAGS = -O3 -Wall -DNDEBUG
LIBPATH = -L. 

#--------------------------------------------------
LIBSOURCES =             \
          CSRmatrix.cc   \
          Graph.cc  

LIBOBJECTS         = $(LIBSOURCES:.cc=.o) 

UTILOBJECTS         = mmio.o mmUtil.o 


#--------------------------------------------------
.SUFFIXES:

%.o : %.cc
	$(CCC) $< -c $(CCFLAGS) $(INCDIRS)

%.o : %.f
	$(F77) $< -c $(CCFLAGS) $(INCDIRS)
#--------------------------------------------------

vpath %.cc $(UTILDIR)

all 		:	lib miniTri
matrix.o	:	CSRmatrix.h 

lib		:	$(LIBOBJECTS) $(UTILOBJECTS)
	ar rvu libSPLA.a $(LIBOBJECTS) $(UTILOBJECTS)

miniTri	:	lib CSRmatrix.h Graph.h
	$(CCC) $(INCDIRS) $(LIBPATH) $(CCFLAGS) -o miniTri.exe miniTri.cc -lSPLA 

clean	:
	rm -f *.o *~ libSPLA.a *.out *.exe
