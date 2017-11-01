#Point this to kokkos installation

# Set the following environment variables:
#    KOKKOS_PATH -- installation directory of Kokkos
#    KK_PATH -- installation directory of Kokkos Kernels

KOKKOS_PATH=/home/mmwolf/software/kokkos/install
KK_PATH=/home/mmwolf/software/kokkos_kernels/install


default: build


include $(KOKKOS_PATH)/Makefile.kokkos

KOKKOS_LIBPATH = $(KOKKOS_PATH)/lib

KK_INC=-I$(KK_PATH)/include


UTILDIRS = ../../../utils/


EXTRA_INC=-I$(UTILDIRS)


LINK = ${CXX}
LINKFLAGS =  
EXE  = triCount


DEPFLAGS = -M

#--------------------------------------------------
LIBSOURCES =             \
          Graph.cpp

LIBOBJECTS         = $(LIBSOURCES:.cpp=.o)
UTILOBJECTS        = mmio.o mmUtil.o 


LIB = -L. -lminiTriGraph
#--------------------------------------------------
vpath %.cc $(UTILDIRS)
vpath %.cpp $(UTILDIRS)
include $(KOKKOS_PATH)/Makefile.kokkos

build: lib $(EXE)

lib		:	$(LIBOBJECTS) $(UTILOBJECTS)
	ar rvu libminiTriGraph.a $(LIBOBJECTS) $(UTILOBJECTS)


$(EXE): lib triCount.o
	$(LINK) -fopenmp $(LINKFLAGS) $(EXTRA_PATH) triCount.o $(LIB) -L$(KOKKOS_LIBPATH) $(KOKKOS_LIBS) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host *.a $(EXE)

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) $(KK_INC) -c $< -o $(notdir $@)
