HDF5INCL= -I/home/leporif/Nbody/Lib/hdf5-1.10.1_lib/include
HDF5LIBS= -L/home/leporif/Nbody/Lib/hdf5-1.10.1_lib/lib
FFTWINCL= -I/home/leporif/Nbody/Lib/fftw-3.3.4_lib/include
FFTWLIBS= -L/home/leporif/Nbody/Lib/fftw-3.3.4_lib/lib
GSLINCL= -I/home/leporif/Nbody/Lib/gsl-1.16_lib/include
GSLLIBS= -L/home/leporif/Nbody/Lib/gsl-1.16_lib/lib
# programming environment
COMPILER     := mpic++
INCLUDE      := -I/home/leporif/Nbody/Lib/LATfield2 $(HDF5INCL) $(FFTWINCL) $(GSLINCL)
LIB          := $(HDF5LIBS) $(FFTWLIBS) $(GSLLIBS) -lfftw3 -lm -lhdf5 -lgsl -lgslcblas 
#-lfftw3 -lm -lhdf5 -lgsl -lgslcblas

# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
#DGEVOLUTION  += -DCHECK_B
#DGEVOLUTION  += -DHAVE_CLASS # requires OPT -fopenmp and LIB -lclass

# further compiler options
OPT          := -O3 -std=c++11 

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

