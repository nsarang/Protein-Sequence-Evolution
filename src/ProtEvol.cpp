//
//  ProtEvol.cpp
//  Protein Sequence Evolution
//
//  Created by Nima Sarang on 2018-02-18.
//  Copyright © 2018 Nima Sarang. All rights reserved.
//

// INCLUDES
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <cassert>
#include <array>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <sstream>
#include <set>
#include <iomanip>
#include "protein.h"
#include "constants.h"
#include "cpplocate.h"
#include "statistics.h"
#include "utility.h"
#include "mapping.h"
#include "db_process.h"
#include "alignment.h"
#include "dDFIRE2.h"
#include "evaluation.h"
#include "init.h"
#include "ai.h"
// #include "helper.cpp"
#define _USE_MATH_DEFINES



// FUNCTIONS



int main(int argc, const char * argv[]) {
    // clocate::SetCurrentWorkingDirectory(clocate::GetDirectoryPath(clocate::GetDirectoryPath(clocate::GetExecutablePath())));
    std::srand(std::time(NULL));


    assert(argc > 0);
    std::string target_path = argv[1];
   // std::cerr << GetSec(target_path) << "\n";
   // std::cerr << File_md5(target_path) << "\n\n";
    // SuperAlgnmntProcess(2000);
    Initialize(target_path);
    // std::cerr << __target_sequence << std::endl;
    auto AI_Obj = AI();

    // std::string rett = AI_Obj.AntColonyOptimization(10, 50, 1, 3, 1, 0.5);
  //  std::string orig = GetSec(target_path);
 //   std::cerr << "ORIG: " << orig << "\n";
  //  std::cout << AI_Obj.Initial_State(orig.size()) << "\n";
    // EVAL(rett, orig);
    return 0;
    // AI_Obj.ParticleSwarmOptimization(500, 10000);

    //std::string init = AI_Obj.Initial_State(), s1, s2;
    // std::cerr << "$$ 1st run $$\n";
    // s1 = AI_Obj.Stochastic_Hill_Climbing(init);
    // s2 = AI_Obj.Simulated_Annealing(init, 1000, 0.99);
    // std::cerr << "HC: " << s1 << "\n" << O_Fitna(s1) << "\n\n" << "SA: " << s2 << "\n" << O_Fitna(s2) << "\n\n";

    // std::cerr << "$$ 2nd run $$\n";
    // s1 = AI_Obj.Stochastic_Hill_Climbing(init);
    // s2 = AI_Obj.Simulated_Annealing(init, 1000, 0.95);
    // std::cerr << "HC: " << s1 << "\n" << O_Fitna(s1) << "\n\n" << "SA: " << s2 << "\n" << O_Fitna(s2) << "\n\n\n";
    /*
    std::cerr << "$$ 3rd run $$\n";
    s1 = AI_Obj.Stochastic_Hill_Climbing(init, 1000);
    s2 = AI_Obj.Simulated_Annealing(init, 10000, 0.95, 1000);
    std::cerr << "HC: " << s1 << "\n" << O_Fitna(s1) << "\n\n" << "SA: " << s2 << "\n" << O_Fitna(s2) << "\n\n\n";
    EVAL(__target_sequence, s1);
    EVAL(__target_sequence, s2);

    std::cerr << "$$ 4th run $$\n";
    s1 = AI_Obj.Stochastic_Hill_Climbing(init, 1000);
    s2 = AI_Obj.Simulated_Annealing(init, 50000, 0.999, 1000);
    std::cerr << "HC:\n" << s1 << "\n" << O_Fitna(s1) << "\n\n" << "SA: " << s2 << "\n" << O_Fitna(s2) << "\n\n\n";
    EVAL(__target_sequence, s1);
    EVAL(__target_sequence, s2);
    */
    // std::cerr << "$$ 5th run $$\n";
    // s1 = AI_Obj.Stochastic_Hill_Climbing(init);
    // s2 = AI_Obj.Simulated_Annealing(init, 100000, 0.99);
    // std::cerr << "HC:\n" << s1 << "\n" << O_Fitna(s1) << "\n\n" << "SA: " << s2 << "\n" << O_Fitna(s2) << "\n\n\n";

    // ScoreEvaluation(PotScore);
    // std::string sl = R"(KEEAAEMAAEMMEAAAAAAEMAAEAMAAARMARAAFAAAASTSELALAAAAAAMWLAAAAAATAEAAAAAALAAAAEAAAAAAAMEMAAAAAAAEAAMAAAAAAMAAMAAAEELAAAAAAAKAAWAAAASAAAAAAMAVAALAEASEAAAESMAAAAALAALAA)";
    // std::cerr << PotScore(sl) << "\n";

    // std::string s2 = R"(KEEAAEMAAEMMEAAAAAAEMAASAMAAARMARAAMAAAAETSELALAAAAAAMALAAAAAATAEAAAAAALAAAAEAAAAAAAMEMAAAAAAAEAAMAAAAAAMAAMAAAEELAAAAAAAKAAWNAAASAAAAAAMAAAALAEAAEAAAEAMAAAAALAALAA)";
    // std::cerr << PotScore(s2) << "\n";

    // std::string s3 = R"(EDEMAEAAAEAEDAAAAAAEAAAEAMAAARAARAAMAALASTESALAAAMAMAAMAAAAMAATAEAAMLAEAAAAAEAAAAAMAADAAAAAAAAEAAAAAAAMMAAAAAAAEAALMLAAAARAAWMAAAAAAAAAAAAAAAAAAAAEMAAEAALLAAAAAALAA)";
    // std::cerr << PotScore(s3) << "\n";
    return 0;
}