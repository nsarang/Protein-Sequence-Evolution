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


// void write_seq_to_pdb(source, dest, seq) {
//     std::ifstream inFile(source);
//     std::ofstream outFile(dest);
//     std::string line;

//     while (getline(inFile, line)) {
//         if (line.size() < 47 || line.substr(0, 4) != "ATOM")
//             continue;

//         std::string aa_name;
//         for (auto&[name, sym] : resName_to_sym) {
//         	if (sym == seq[aa_pos])
//         }

//         last = line.substr(22, 5);
//         retVec.emplace_back(AminoAcid(
//                                 resName_to_sym.count(line.substr(17, 3)) ?  line.substr(17, 3) : "UNK",
//                                 ' ',
//                                 std::stod(utility::trim(line.substr(30, 8))),
//                                 std::stod(utility::trim(line.substr(38, 8))),
//                                 std::stod(utility::trim(line.substr(46, 8)))
//                             ));

//         retVec.back().symbol = resName_to_sym[retVec.back().name];
//     }
// }


int main(int argc, const char * argv[]) {
	// Parameters
	assert(argc == 3);
	std::string pdb_path = argv[1],
				profile_dir = argv[2];


	// Load profiles
	auto prot = Protein(pdb_path, 1 + 2 + 4 + 8);
	auto profile = ProteinProfile(prot);
	
	profile.Read_FromFile(profile_dir, 1 + 2 + 4 + 8);
	if (profile.RemainingProfiles() != 0)
	{
		auto vecDB = utility::CATH_ListFiles(db_CATH);
		profile.Find_Homologous_Proteins( vecDB, true );
		profile.CalculateProfiles( profile.RemainingProfiles(), true );
		profile.Write_ToFile(profile_dir, 1 + 2 + 4 + 8 + 32);
	}


	// Output scores
	std::cout << "OUTPUT: ";
	DeepAI DI;
	DI.GenerateData(profile, prot, std::cout, " ");
}