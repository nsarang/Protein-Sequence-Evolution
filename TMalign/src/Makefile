FC=gfortran-7
FFLAGS= -O3 -ffast-math -lm
CC=g++-7
CFLAGS=-static

JC=javac
JAR=jar

FORTRAN_PROGRAMS=TMalign

all: ${FORTRAN_PROGRAMS}

TMalign: TMalign.f
	${FC} ${FFLAGS} -o ../TMalign TMalign.f


install: all
	strip ${FORTRAN_PROGRAMS}

clean: 
	rm -f ${FORTRAN_PROGRAMS} 
