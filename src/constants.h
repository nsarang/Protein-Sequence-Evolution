#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <string>


// PATHS
const std::string db_CATH = "./CATH/";
const std::string db_Profiles = "./Profiles/";
const std::string db_Align = "./Alignments/";
const std::string ex_STRIDE = "./STRIDE/stride";
const std::string ex_TMALIGN = "./TMalign/TMalign";
const std::string sFileStartSym = "~^&%@";

// MATH
#define M_EPS 1e-15


// PROTEIN
const int MAX_PROT_LEN = 2000,
          AA_TYPE = 20;


// POT STATISTIC
const int Pot_S_Constant = 2;


// ALIGNMENT REPOSITORY
const double DIST_CUTOFF = .5,
             GAP_PENALTY = 20,
             GAP_NUMBER = 99,
             ALIGN_SC_CUTTOFF = .4;
const int FRAG_MIN_LEN = 4;


// WEIGHTS
const int WGT_NUM = 6;


// dDFIRE
const int MAX_TYPES = 200,
          MAX_BIN = 40;
const double dDFIRE_COEF = -4.11823;



#endif // __CONSTANTS_H__