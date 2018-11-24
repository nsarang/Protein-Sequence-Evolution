#include <cmath>
#include <array>
#include <vector>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <functional>
#include <dirent.h>
#include <unistd.h>
#include <tuple>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "statistics.h"
#include "protein.h"
#include "protein_profile.h"
#include "evaluator.h"


int main(int argc, const char * argv[]) {
    std::srand(std::time(NULL));

    assert(argc > 1);
    std::string target_path = argv[1];

    auto target = Protein(target_path, 1 + 2 + 4 + 8);
 //  for (int i = 0; i < 20; ++i) {
 //      std::cerr <<target.aAA_Freqs[i] << " \n"[i == target.length()-1]; 
 //  }
  // return 0;
//        for (int i = 0; i < target.length(); ++i) {
//        std::cerr << target.aSecondary_Structure[i] << " \n"[i == target.length()-1]; 
//    }
    std::cerr << target.Get_Sequence() << "\n";

    auto profile = ProteinProfile(target);

    std::vector<Protein> pt;
    auto vecDB = CATH_ListFiles(db_CATH);
 //   for (auto path : vecDB)
   //     system(("md5 -r " + path).c_str());
  //      pt.push_back(Protein(path));


    profile.Find_Homologous_Proteins(vecDB, DIST_CUTOFF, true);
    profile.CalculateProfiles(true, true, true, true);
    profile.Write_ToFile(true);

}