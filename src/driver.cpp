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
    std::cerr << utility::CountFilesInDir(".") << "\n";

    auto target = Protein(target_path, 1 + 2 + 4 + 8);
    auto profile = ProteinProfile(target);

 //  for (int i = 0; i < 20; ++i) {
 //      std::cerr <<target.aAA_Freqs[i] << " \n"[i == target.length()-1]; 
 //  }
  // return 0;
//        for (int i = 0; i < target.length(); ++i) {
//        std::cerr << target.aSecondary_Structure[i] << " \n"[i == target.length()-1]; 
//    }
    std::cerr << target.Get_Sequence() << "\n";


 //   for (auto path : vecDB)
   //     system(("md5 -r " + path).c_str());
  //      pt.push_back(Protein(path));
    auto vecDB = utility::CATH_ListFiles(db_CATH);


    profile.Find_Homologous_Proteins(vecDB, DIST_CUTOFF, true);
    profile.CalculateProfiles(1 + 2 + 4 + 8);
    profile.Write_ToFile(true);

}