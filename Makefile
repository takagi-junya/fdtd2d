FC		= ftn
CC		= cc

FFLAGS		= -ipo -O3 -no-prec-div -g -static -fopenmp -fp-model fast=2 -xMIC-AVX512 -fpp -Dssurf=3
CFLAGS		= -ipo -O3 -no-prec-div -g -static -fopenmp -fp-model fast=2 -xMIC-AVX512
LFLAGS		= -ipo -O3 -no-prec-div -g -static -fopenmp -fp-model fast=2 -xMIC-AVX512 -fpp -Dssurf=3

LDFLAGS		= 
HDF5DIR		= ./
FFTWDIR		= ./
MPIDIR = ./

LINKER		= $(FC)

MODS 	= constants.o tfsf_gausian.o tfsf_dgausian.o tfsf_sin.o tfsf_RCP.o tfsf_LCP.o pml2d.o hdfio.o output.o farfield.o directivity.o

FOBJS	= fdtd2d.o setup.o PEC.o dielec.o func.o efield.o hfield.o ADE.o ADEg.o JEC.o JECg.o EOM.o EOMg.o current.o currentEOM.o velocity.o comm2d.o finalize.o

OBJS	= $(MODS) $(FOBJS) 

.SUFFIXES: .F90

fdtd:	$(OBJS)
		$(LINKER) $(OBJS) $(LFLAGS) $(LDFLAGS) -o $@

.F90.o:
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -I$(MPIDIR) -c $< -o $@
   
$(MODS):%.o:	%.F90 
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -c $< -o $@
$(KMODS):%.o:	%.F90
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -c $< -o $@
clean:;
		rm *.o; rm *.mod
