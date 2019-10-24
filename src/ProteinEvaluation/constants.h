#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <string>


// PATHS
const std::string db_CATH = "./data/CATH/";
const std::string db_Profiles = "./data/Profiles/";
const std::string db_DFIRE2 = "./data/dDFIRE.lib";
const std::string fl_CATH_BL = "./data/cath_blacklist.txt";
const std::string ex_STRIDE = "./bin/stride";
const std::string ex_TMALIGN = "./bin/TMalign";
const std::string sFileStartSym = "~^&%@";

// MATH
#define M_EPS 1e-15


// THREAD
const int FALLBACK_THRD_NUM = 8;


// PROTEIN
//const int MAX_PROT_LEN = 2000,
//          AA_TYPE = 20;


// POT STATISTIC
const int Pot_S_Constant = 2;


// ALIGNMENT REPOSITORY
const double DIST_CUTOFF = .8,
             GAP_PENALTY = 15,
             GAP_NUMBER = 99,
             ALIGN_SC_CUTTOFF = .5;
const int MAX_FAM_SIZE = 20,
		  FRAG_MIN_LEN = 4;


// WEIGHTS
// const int WGT_NUM = 6;


// dDFIRE
// const int MAX_TYPES = 200,
//          MAX_BIN = 40;
const double dDFIRE_COEF = -4.11823;



#endif // __CONSTANTS_H__