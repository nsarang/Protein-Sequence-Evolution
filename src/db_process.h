#ifndef __DB_PROCESS_H__
#define __DB_PROCESS_H__


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
#include <dirent.h>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "statistics.h"



// PATHS
extern const std::string db_CATH;


// DATABASE FILE LIST
extern std::vector<std::string> db_protein_list;


// SCORE PROFILES
extern std::array<std::array<double, 7>, 20> burial_potential,
    structure_potential;


// POT STATISTICS
extern std::array<double, 20> pot_bar, pot_stdev;
extern const int Pot_S_Constant;


// AA FREQUENCY
extern std::array<long long, 20> AA_frequency;
extern std::array<double, 20> AA_freq_mean, AA_freq_sd;





// FUNCTIONS
void Burial_PotStatistic(std::string fPath);

void SecondaryStructure(std::string fPath);

void BurSecPot_Eval();

void DB_ListFiles(std::string path = db_CATH);

void PreprocessThread(int s, int t);

void DB_PreprocMultiThrd(int thread_num = 12);


#endif // __DB_PROCESS_H__
