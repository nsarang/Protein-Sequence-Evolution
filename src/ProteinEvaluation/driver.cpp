#include <iostream>
#include <algorithm>
#include <random>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "statistics.h"
#include "protein.h"
#include "protein_profile.h"
#include "evaluator.h"
#include "deep.h"
using namespace std;


int main(int argc, const char * argv[]) {
	// assert(argc > 1);
	// std::string target_path = argv[1];
	// auto target = Protein(target_path, 1 + 2 + 4 + 8);


	auto vecDB = utility::CATH_ListFiles(db_CATH);
	DeepAI DI;
	// DI.CheckFamilies("./FamilyProfiles/");
	DI.GenerateDataset("./FamilyProfiles/", "dataset.csv", 10);


	// std::ofstream trFile("ole.csv");


	// for (int i = 0; i < 1000; ++i) {
	// 		auto pr = ProteinProfile(Protein("./CATH/1y6dA00"));
	// 		pr.Read_FromFile(db_Profiles, 1 + 2 + 4 + 8);
	// 	auto po = Protein(vecDB[10]);
	// 	DI.GenerateData(pr, po, trFile);
	// 	std::cerr << i << "\n";
	// }

	// DI.GenerateDataset("./FamilyProfiles/", "dataset.csv", 10);
	// auto vecProfiles = DI.PrepareProfiles("./FamilyProfiles/");
	//DI.GenerateDataset(vecProfiles, "dataset.csv", 10);


/*
	std::string fPath1 = "./CATH/2nvnA00", fPath2 = "./CATH/1dznA01";
	auto prot = Protein(fPath1);
	auto profile = ProteinProfile(prot);
	profile.Read_FromFile();
	auto prot2 = Protein(fPath2, 1 + 2 + 4 + 8);
	std::ofstream o("ole.txt");
	DI.GenerateScores(profile, prot2, o);
	*/

	/*
	decltype(vecDB) vecSample;
	std::sample(vecDB.begin(), vecDB.end(),
	            std::back_inserter(vecSample),
	            1000, std::mt19937{std::random_device{}()});
	*/

	/*
	for (auto fPath : vecDB) {
		std::cerr << fPath << "\n";
		auto prot = Protein(fPath);
		auto profile = ProteinProfile(prot);

		profile.Read_FromFile();
		if (profile.CalculateRemainingProfiles(vecDB, true))
			profile.Write_ToFile();
	}
	*/
}