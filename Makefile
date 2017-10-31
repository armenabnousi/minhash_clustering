MRMPI_LIB := /home/armen.abnousi/libs/mrmpi-7Apr14

CC =            mpicc
CPP =           mpic++
CPPFLAGS =      -g -O -ggdb -std=c++14 -I$(MRMPI_LIB)/src
LINK =          mpic++
LINKFLAGS =     -g -O -Wall
USRLIB =        $(MRMPI_LIB)/src/libmrmpi_mpicc.a
SYSLIB =        -lpthread


export MRMPI_LIB
export CC
export CPP
export CPPFLAGS
export LINK
export LINKFLAGS
export USRLIB
export SYSLIB

all:
	cd src && ${MAKE}
	cd graph_converter && ${MAKE}
	cd fvalue_evaluator && ${MAKE}

clean:
	rm -f */*.o */*/*.o
