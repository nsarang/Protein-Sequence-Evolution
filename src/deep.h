#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include "protein.h"
#include "protein_profile.h"
#include "evaluator.h"
#include "constants.h"



class DeepAI {
public:
	std::vector<ProteinProfile> PrepareProfiles(std::string sFDir);

	void GenerateDataset(std::vector<ProteinProfile>& vecProfiles,
	                     std::string fCSV,
	                     int maxPerProt);

	void GenerateScores(ProteinProfile& profile, Protein& prot,
	                    std::ofstream& outFile, std::string sep = ",");

	double PairAlgnScore(std::string PDB_1, std::string PDB_2);
};

