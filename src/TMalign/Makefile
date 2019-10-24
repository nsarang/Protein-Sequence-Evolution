FC=gfortran-8
FFLAGS= -O3 -ffast-math -lm
CC=g++-8
CFLAGS=-static

JC=javac
JAR=jar

FORTRAN_PROGRAMS=TMalign

all: ${FORTRAN_PROGRAMS}

TMalign: TMalign.f
	${FC} ${FFLAGS} -o "$(BINDIR)/TMalign" TMalign.f


install: all
	strip ${FORTRAN_PROGRAMS}

clean: 
	rm -f ${FORTRAN_PROGRAMS} 
