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
	
	/*
	decltype(vecDB) vecSample;
	std::sample(vecDB.begin(), vecDB.end(),
	            std::back_inserter(vecDB1),
	            1000, std::mt19937{std::random_device{}()});
	*/

	for (auto fPath : vecDB) {
		std::cerr << fPath << "\n";
		auto prot = Protein(fPath);
		auto profile = ProteinProfile(prot);

		profile.Read_FromFile(db_Profiles, 1 + 2 + 4 + 8);
		if (profile.RemainingProfiles() != 0)
		{
			profile.Read_FromFile(db_Profiles, 32);
			if (profile.FamilySize() == 0) {
				profile.Find_Homologous_Proteins(vecDB);
		        }
 			profile.CalculateProfiles( profile.RemainingProfiles() );
			profile.Write_ToFile(db_Profiles, 1 + 2 + 4 + 8 + 32);
		}
	}
}
