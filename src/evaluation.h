#ifndef __EVALUATION_H__
#define __EVALUATION_H__


#include <iostream>
#include "mapping.h"
#include "constants.h"
#include "db_process.h"
#include "alignment.h"
#include "dDFIRE2.h"
#include "init.h"



// FUNCTIONS
double PotScore(std::string& target_sequence);

double BurialScore(std::string& target_sequence);

double SecondaryStructScore(std::string& target_sequence);

double AlignmentScore(std::string& target_sequence);

double FrequencyScore(std::string& target_sequence);

double O_Fitna(std::string target_sequence, double ret[] = NULL);

double eNormalaize(double e, double lB, double uB);



#endif // __EVALUATION_H__
