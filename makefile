# programming environment
COMPILER     := mpic++
INCLUDE      := -I/home/jelicciz/LATfield2# $(HDF5INCL) $(FFTWINCL) $(GSLINCL)
LIB          := -lfftw3 -lm -lhdf5 -lgsl -lgslcblas 
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

