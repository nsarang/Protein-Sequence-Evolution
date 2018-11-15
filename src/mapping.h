#ifndef __MAPPING_H__
#define __MAPPING_H__

#include <unordered_map>
#include <array>
#include <string>

// MAPPINGS
extern std::unordered_map<std::string, char> residue_map;

extern std::unordered_map<char, int> symbol_map;

extern char index_to_symbol[];



// CLASSIFICATIONS
extern std::array<int, 6> burial_classification;

extern std::unordered_map<char, int> sec_structure_classification;

extern std::unordered_map<char, int> AA_classification;


#endif // __MAPPING_H__
