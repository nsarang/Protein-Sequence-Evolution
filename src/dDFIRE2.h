#ifndef _DDFIRE2_H
#define _DDFIRE2_H

#include <string>
#include <array>
#include <unordered_map>
#include <fstream>
#include "constants.h"
#include "mapping.h"
#include "utility.h" 


extern std::unordered_map<std::string, int> atom_map;
extern double edDFIRE[MAX_TYPES][MAX_TYPES][MAX_BIN];
extern int num_types, num_bins;



// FUNCTIONS
void dDFIRE_ReadLib(std::string libPath);

double dDFIRE_CFE(std::string& target_sequence,
                  std::array<std::array<double, MAX_PROT_LEN>, MAX_PROT_LEN> &target_atom_distance);


#endif // _DDFIRE2_H
