#Inverse Protein Folding

CC= gcc-8
CXX=g++-8
CPPFLAGS = -std=c++17 -c -Wall -pedantic -O3
LDLIBS = -lm -lgsl -lgslcblas -L. -pthread
RM=rm -f

SRCS=driver.cpp protein.cpp protein_profile.cpp evaluator.cpp utility.cpp mapping.cpp dDFIRE2.cpp statistics.cpp
OBJS=$(subst .cpp,.o,$(SRCS))



all: $(OBJS)
	$(CXX) -o ../NewProtEvol $(OBJS) $(LDLIBS)


depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend