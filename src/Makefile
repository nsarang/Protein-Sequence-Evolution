#Protein Sequence Evolution

PROG=NewProtEvol
CC=gcc-8
CXX=g++-8
CXXFLAGS+=-std=c++17 -g -Wall -pedantic
#CXXFLAGS+= `gsl-config --cflags`
LDFLAGS+= -lstdc++fs -pthread
#LDFLAGS+= `gsl-config --libs` 
RM=rm -f

SRCS=driver.cpp protein.cpp protein_profile.cpp evaluator.cpp utility.cpp mapping.cpp dDFIRE2.cpp statistics.cpp deep.cpp
OBJS=$(subst .cpp,.o,$(SRCS))



all: $(OBJS)
	$(CXX) -o ../$(PROG) $(OBJS) $(LDFLAGS)


depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend