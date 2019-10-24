/*

THIS FILE IS NOT COMPATIBLE WITH THE NEWER VERSION

*/




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
#define _USE_MATH_DEFINES



// FUNCTIONS
void Generate_Random() {
    double ret[6];
    std::srand(time(NULL));
    int c = 50000;

    for (int j = 0; j < 1; ++j) {
        std::string vmax, s, ss;
        double vm = 0;
        for (int k = 0; k < PROT_LEN; ++k)
            s += index_to_symbol[rand() % 20];
        for (int i = 0; i < c; ++i) {
            ss = s;
            int cnt = rand() % 7;
            while (cnt--) {
                int pp = rand() % PROT_LEN;
                ss[pp] = index_to_symbol[rand() % 20];
            }

            // for (int i = 0; i < 5; ++i)
            //     seq[i].push_back(ret[i]);
            // seq[5].push_back(zz);
            std::cerr << "2\n";
            double zz = O_Fitna(ss, ret);
            std::cerr << "3\n";

            std::cerr << "OLE" << "\n";
            if (zz >= vm) {
                vm = zz;
                vmax = ss;
                s = ss;
            }

            std::cerr << zz << "\n";
            std::cerr << ss << "\n\n";
        }
        std::cerr << vm << " " << vmax << "\n";
    }
}


template <typename Func>
void ScoreEvaluation(Func tFitness) {
    std::vector<double> vecSC, vecLen;
    DB_ListFiles();

    double m, sd;
    double mxx = -10000, mnn = 10000;
    std::string fX, fN, seqX, seqN;

    for (auto &fName : db_protein_list)
    {
        __target_sequence = "";
        TargetProteinProfile(db_CATH + fName);
        double sc = tFitness(__target_sequence);
        vecSC.push_back(sc);

        if (mxx < sc) {
            mxx = sc;
            fX = fName;
            seqX = __target_sequence;
        }
        if (mnn > sc) {
            mnn = sc;
            fN = fName;
            seqN = __target_sequence;
        }
    }
    vecSC = OutlierElimination_IQR(vecSC);

    Mean_SD(vecSC, m, sd);
    std::cerr << m << " " << sd << "\t" << mxx << "  " << mnn << "\n";
    std::cerr << fX << "\n" << seqX << "\n\n" << fN << "\n" << seqN << "\n";

    std::sort(vecSC.begin(), vecSC.end());
    std::ofstream outFile("FitnessEval.txt");// outFile1("Stat_dDFire2.txt");
    for (auto d : vecSC)
        outFile << d << "\n";
    // for (auto d : vecLen)
    //     outFile1 << d << "\n";
}


void FrequencyEvaluation() {
    std::vector<double> vecSC[20];
    DB_ListFiles();

    for (auto &fName : db_protein_list)
    {
        std::array<double, 20> freq{};
        __target_sequence = "";

        TargetProteinProfile(db_CATH + fName);

        for (auto c : __target_sequence)
        {
            if (symbol_map.count(c) > 0)
            {
                freq[symbol_map[c]]++;
            }
        }
        // std::cerr << "A\n";
        for (int i = 0; i < 20; ++i)
        {
            freq[i] /= PROT_LEN;
            double sc = std::abs((freq[i] - AA_freq_mean[i]) / AA_freq_sd[i]);
            vecSC[i].push_back(sc);
        }
    }
    std::cerr << "D\n";
    std::ofstream wFile(db_Scores + "weights.txt", std::ios_base::app);
    for (int i = 0; i < 20; ++i)
    {
        std::sort(vecSC[i].begin(), vecSC[i].end());
        vecSC[i] = OutlierElimination_IQR(vecSC[i], 1.5);
        std::ofstream outRes(std::string() + "AA_SC-" + index_to_symbol[i] + ".txt");
        for (auto s : vecSC[i])
            outRes << s << "\n";
        wFile << "WGHT" << index_to_symbol[i] << " " << std::fixed << std::setprecision(4) <<
              vecSC[i].back() << " " << vecSC[i].front() << "\n";
    }
}



void dDFIRE_Regression() {
    std::vector<double> vecSC, vecLen;
    std::srand(time(NULL));
    DB_ListFiles();

    for (auto &fName : db_protein_list)
    {
        __target_sequence = "";
        TargetProteinProfile(db_CATH + fName);
        if (fName == "1gjjA00")
            continue;
        double sc = dDFIRE_CFE(__target_sequence, target_atom_distance);
        vecSC.push_back(sc);
        vecLen.push_back(__target_sequence.length());
    }
}


void SuperAlgnmntProcess(int sz) {
    DB_ListFiles();

    std::set<int> idxs;
    while(idxs.size() != sz)
        idxs.insert(rand() % db_protein_list.size());

    std::ofstream outFile("ProtList.txt");
    for (auto i : idxs)
        outFile << db_protein_list[i] << "\n";
    outFile.close();

    for (auto i : idxs)
        GetAlignments(db_protein_list[i], false);
}

std::string GetSec(std::string target) {
    std::vector<AminoAcid> vec_aa;
    ParsePDB(target, vec_aa);
    // std::cerr << target << " " << vec_aa.size() << "\n";
    std::string sec;
    for (auto aa : vec_aa)
        sec += residue_map[aa.name];//, std::cerr << aa.name << "\n";
    return sec;
}


void EVAL(std::string target, std::string pred) {
    int match = 0;
    for (int i = 0; i < target.length(); ++i)
        match += AA_classification[target[i]] == AA_classification[pred[i]];
    std::cerr << "AQ: " << (double)(match * 100) / target.length() << "\n";
}

/*
int main(int argc, const char * argv[]) {
    // clocate::SetCurrentWorkingDirectory(clocate::GetDirectoryPath(clocate::GetDirectoryPath(clocate::GetExecutablePath())));
    std::srand(std::time(NULL));


    assert(argc > 0);
    std::string target_path = argv[1];
   // std::cerr << GetSec(target_path) << "\n";
   // std::cerr << File_md5(target_path) << "\n\n";
    // SuperAlgnmntProcess(2000);
    // Initialize(target_path);
    // std::cerr << __target_sequence << std::endl;
    // auto AI_Obj = AI();


    // std::string rett = AI_Obj.AntColonyOptimization(10, 50, 1, 3, 1, 0.5);
    std::string orig = GetSec(target_path);
    std::cerr << "ORIG: " << orig << "\n";
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
//    return 0;
//}
