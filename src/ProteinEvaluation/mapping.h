#ifndef __MAPPING_H__
#define __MAPPING_H__

#include <unordered_map>
#include <array>
#include <string>



// MAPPINGS
extern std::unordered_map<std::string, char> resName_to_sym;

extern std::unordered_map<char, int> sym_to_idx;

extern char idx_to_sym[];



// CLASSIFICATIONS
extern std::array<int, 6> solvent_classes;

extern std::unordered_map<char, int> sec_classes;

extern std::unordered_map<char, int> algn_classes;


#endif // __MAPPING_H__
