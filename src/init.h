#ifndef __INIT_H__
#define __INIT_H__


#include <array>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "db_process.h"
#include "alignment.h"
#include "dDFIRE2.h"


// TARGET PROTEIN PROFILE
extern int PROT_LEN;
extern std::array<int, MAX_PROT_LEN> target_secondary_structure,
    target_solvent_accessibility;
extern std::array<std::array<double, MAX_PROT_LEN>, MAX_PROT_LEN> target_atom_distance;
extern std::string __target_sequence;
extern std::array<double, 20> target_AA_frequency;


// WEIGHTS
extern std::array<std::array<double, 3>, WGT_NUM> sc_weight;
extern std::array<std::array<double, 2>, 20> freq_bounds;



// FUNCTIONS
void Initialize(std::string target_path);

void TargetProteinProfile(std::string fName);

void GetScores();

void GetWeights();

void GetAlignments(std::string target, bool prcss_frgmnts = true);

#endif // __INIT_H__
