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

#Set environment variable HPX_LOCATION to the location of the HPX installation
export PKG_CONFIG_PATH := $(PKG_CONFIG_PATH):$(HPX_LOCATION)/lib/pkgconfig

HPX_CFLAGS := `pkg-config --cflags hpx_application`
HPX_LINKERFLAGS := `pkg-config --cflags --libs hpx_application`

UTILDIR = ../../utils/
INCDIRS = -I. -I$(UTILDIR)
CCC = g++                                                                                                                                   
CCFLAGS = -O3 -Wall -DNDEBUG $(HPX_CFLAGS)
LINKERFLAGS = -O3 -Wall -DNDEBUG $(HPX_LINKERFLAGS)
LIBPATH = -L.

#--------------------------------------------------
LIBSOURCES =             \
          CSRMatrix.cc   \
          Graph.cc  

LIBOBJECTS         = $(LIBSOURCES:.cc=.o) 
UTILOBJECTS	   = mmio.o mmUtil.o
#--------------------------------------------------

.SUFFIXES:

%.o : %.cc
	$(CCC) $< -c $(CCFLAGS) $(INCDIRS)

%.o : %.f
	$(F77) $< -c $(CCFLAGS) $(INCDIRS)


#--------------------------------------------------

vpath %.cc $(UTILDIR)

all 		:	lib triangleEnum
matrix.o	:	CSRMatrix.h 

lib		:	 $(LIBOBJECTS) $(UTILOBJECTS)
	ar rvu libSPLA.a $(LIBOBJECTS) $(UTILOBJECTS)

triangleEnum	:	lib CSRMatrix.h Graph.h
	$(CCC) $(INCDIRS) $(LIBPATH) -o miniTri.exe miniTri.cc -lSPLA $(LINKERFLAGS)


clean	:
	rm -f *.o *~ libSPLA.a *.out miniTri.exe
