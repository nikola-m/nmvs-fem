#
# Makefile for fem2d_elastica_cg program
#         

F90FLAGS = -O2 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      precision.f90 \
      matrix.f90 \
      quadrature.f90 \
      shape_functions.f90 \
      fem_module.f90 \
      fem_mesh.f90 \
      output.f90 \
      fem2d_elastica_cg.f90



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: fem2d_elastica_cg

fem2d_elastica_cg: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o fem2d_elastica_cg

.PHONY: clean
clean:
	@rm  *.o fem2d_elastica_cg

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
