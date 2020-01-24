#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <set>
#include "protein.h"
#include "protein_profile.h"
#include "evaluator.h"
#include "constants.h"
#include "ThreadPool.h"
#include "utility.h"



class Dataset {
public:
	void CheckFamilies(std::string sFDir);

	void GenerateCASPDataset(std::string caspTargetsDir,
							 std::string serverPredictionsDir,
							 std::string fCSVOutputPath);

	void GenerateDataset(std::string sFamDir,
						 std::string fCSV,
						 int nFamCutoff);

	void GenerateData(ProteinProfile& profile, Protein prot,
	                  std::ostream& outFile, std::string sep = ",");

	double PairAlgnScore(std::string PDB_1, std::string PDB_2);
};

