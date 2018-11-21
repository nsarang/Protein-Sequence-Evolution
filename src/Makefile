#Inverse Protein Folding

CC= g++-8
CFLAGS = -std=c++17 -c -Wall -pedantic -O3
LDFLAGS = -lm -lgsl -lgslcblas -L.
ALL.O = ProtEvol.o cpplocate.o statistics.o utility.o db_process.o mapping.o alignment.o dDFIRE2.o init.o evaluation.o ai.o protein.o
ALL.H = cpplocate.h statistics.h utility.h mapping.h db_process.h constants.h alignment.h dDFIRE2.h init.h evaluation.h ai.h protein.h
ALL.CPP = ProtEvol.cpp cpplocate.cpp statistics.cpp utility.cpp db_process.cpp mapping.cpp alignment.cpp dDFIRE2.cpp init.cpp evaluation.cpp ai.cpp protein.cpp

all: $(ALL.O)
	$(CC) -o ../ProteinSequenceEvolution $(ALL.O) $(LDFLAGS)

new: protein.o protein_profile.o evaluator.o
	$(CC) -o ../NewEvol protein.o protein_profile.o evaluator.o $(LDFLAGS)

protein.o: protein.cpp protein.h
	$(CC) $(CFLAGS) protein.cpp

protein_profile.o: protein_profile.cpp protein_profile.h
	$(CC) $(CFLAGS) protein_profile.cpp

evaluator.o: evaluator.cpp evaluator.h
	$(CC) $(CFLAGS) evaluator.cpp

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


ai.o: ai.cpp ai.h
	$(CC) $(CFLAGS) ai.cpp

ProtEvol.o: ProtEvol.cpp $(ALL.H)
	$(CC) $(CFLAGS) ProtEvol.cpp


clean: 
	rm -rf *.o ProteinSequenceEvolution
