#include "mapping.h"



// MAPPINGS
std::unordered_map<std::string, char> residue_map = {
    {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'},
    {"ASX", 'B'}, {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'},
    {"GLX", 'Z'}, {"GLY", 'G'}, {"HIS", 'H'}, {"ILE", 'I'},
    {"LEU", 'L'}, {"LYS", 'K'}, {"MET", 'M'}, {"PHE", 'F'},
    {"PRO", 'P'}, {"SER", 'S'}, {"THR", 'T'}, {"TRP", 'W'},
    {"TYR", 'Y'}, {"VAL", 'V'}, {"UNK", 'X'}, {"SEC", 'U'},
    {"PYL", 'O'}
};

std::unordered_map<char, int> symbol_map = {
    {'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4},
    {'Q', 5}, {'E', 6}, {'G', 7}, {'H', 8}, {'I', 9},
    {'L', 10}, {'K', 11}, {'M', 12}, {'F', 13}, {'P', 14},
    {'S', 15}, {'T', 16}, {'W', 17}, {'Y', 18}, {'V', 19}
};

char index_to_symbol[] = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
};



// CLASSIFICATIONS
std::array<int, 6> burial_classification = {27, 34, 40, 47, 55, 66};

std::unordered_map<char, int> sec_structure_classification = {
    {'H', 0}, {'G', 1}, {'I', 2}, {'E', 3},
    {'B', 4}, {'T', 5}, {'C', 6}
};

std::unordered_map<char, int> AA_classification = {
    {'A', 0}, {'V', 0}, {'L', 0}, {'I', 0}, {'M', 0}, {'C', 0},
    {'G', 1}, {'S', 1}, {'T', 1},
    {'D', 2}, {'E', 2},
    {'N', 3}, {'Q', 3},
    {'R', 4}, {'K', 4},
    {'P', 5}, {'F', 5}, {'Y', 5}, {'W', 5},
    {'H', 6}
};
