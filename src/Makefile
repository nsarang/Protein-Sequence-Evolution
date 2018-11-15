#Inverse Protein Folding

CC= g++-7
CFLAGS = -c -Wall -pedantic -O3
LDFLAGS = -lm -lgsl -lgslcblas -L.
ALL.O = ProtEvol.o cpplocate.o statistics.o utility.o db_process.o mapping.o alignment.o dDFIRE2.o init.o evaluation.o ai.o
ALL.H = cpplocate.h statistics.h utility.h mapping.h db_process.h constants.h alignment.h dDFIRE2.h init.h evaluation.h ai.h
ALL.CPP = ProtEvol.cpp cpplocate.cpp statistics.cpp utility.cpp db_process.cpp mapping.cpp alignment.cpp dDFIRE2.cpp init.cpp evaluation.cpp ai.cpp

all: $(ALL.O)
	$(CC) -o ../ProteinSequenceEvolution $(ALL.O) $(LDFLAGS)


cpplocate.o: cpplocate.cpp cpplocate.h
	$(CC) $(CFLAGS) cpplocate.cpp


statistics.o: statistics.cpp statistics.h
	$(CC) $(CFLAGS) statistics.cpp


utility.o: utility.cpp utility.h
	$(CC) $(CFLAGS) utility.cpp


db_process.o: db_process.cpp db_process.h
	$(CC) $(CFLAGS) db_process.cpp


mapping.o: mapping.cpp mapping.h
	$(CC) $(CFLAGS) mapping.cpp


alignment.o: alignment.cpp alignment.h
	$(CC) $(CFLAGS) alignment.cpp


dDFIRE2.o: dDFIRE2.cpp dDFIRE2.h
	$(CC) $(CFLAGS) dDFIRE2.cpp


init.o: init.cpp init.h
	$(CC) $(CFLAGS) init.cpp


evaluation.o: evaluation.cpp evaluation.h
	$(CC) $(CFLAGS) evaluation.cpp


ai.o: ai.cpp ai.h
	$(CC) $(CFLAGS) ai.cpp


ProtEvol.o: ProtEvol.cpp $(ALL.H)
	$(CC) $(CFLAGS) ProtEvol.cpp


clean: 
	rm -rf *.o ProteinSequenceEvolution
