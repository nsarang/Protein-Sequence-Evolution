#ifndef __ALIGNMENT_H__
#define __ALIGNMENT_H__


#include <cmath>
#include <array>
#include <vector>
#include <thread>
#include <mutex>
#include <tuple>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include "mapping.h"
#include "utility.h"
#include "constants.h"
#include "db_process.h"


// ALIGNMENT REPOSITORY
extern std::vector<std::vector<std::vector<std::string>>> matAMR;


// SCORE PROFILES
extern int PROT_LEN;
extern int TEMPLATE_COUNT;
extern std::array<std::array<double, 7>, MAX_PROT_LEN> alignment_potential;



// FUNCTIONS
void StructureAlignment(std::string target, std::string fPath, bool prcss_frgmnts);

void Alignment_Eval(std::string target);

void AlignmentThread(std::string &target, bool prcss_frgmnts, int s, int t);

void DB_AlignMultiThrd(std::string target, bool prcss_frgmnts, int thread_num = 12);


#endif // __ALIGNMENT_H__
