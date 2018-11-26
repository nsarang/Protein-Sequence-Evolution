
#include <iostream>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "statistics.h"
#include "protein.h"
#include "protein_profile.h"
#include "evaluator.h"





int main(int argc, const char * argv[]) {
    assert(argc > 1);
    std::string target_path = argv[1];

    auto target = Protein(target_path, 1 + 2 + 4 + 8);
    std::cerr << target.Get_Sequence() << "\n";

    auto profile = ProteinProfile(target);

    auto vecDB = utility::CATH_ListFiles(db_CATH);
    profile.Find_Homologous_Proteins(vecDB, DIST_CUTOFF, true);
    profile.CalculateProfiles(1 + 2 + 4 + 8, true);
    profile.Write_ToFile(true);
}