#ifndef __UTILITY_H__
#define __UTILITY_H__


#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <array>
#include <algorithm>
#include <string>
#include <dirent.h>
#include <cerrno>
#include <cstring>
#include "cpplocate.h"
#include "mapping.h"
#include "subprocess.hpp"


#define RESET   "\033[0m"
#define BOLDRED     "\033[1m\033[31m"
#define BOLDGREEN   "\033[1m\033[32m"
#define BUFFERSIZE 1024




// TYPE DEFINITIONS
/*
struct AminoAcid {
    AminoAcid(std::string iname, double x = 0, double y = 0, double z = 0,
              int count = 0)
        : name( iname ), neighbour_count( count ), cords( std::tie(x, y, z)) {}

    std::string name;
    int neighbour_count;
    std::tuple<double, double, double> cords;
};
*/


// FUNCTIONS
std::vector<std::string> CATH_ListFiles(std::string path);

double dist(std::tuple<double, double, double> &t1, std::tuple<double, double, double> &t2);

bool IsStandardAA(std::string abrv);

void Progress_Indicator(std::string text, long long current, long long total);

int system_call_err(std::string command, std::string& stdout);

std::string system_call(std::string command);

void ltrim(std::string &s);

void rtrim(std::string &s);

std::string trim(std::string s);

bool FileExists(std::string fName);

size_t FileSize(std::istream &isObj);

std::string FileBasename(std::string filename);

std::string File_md5(std::string fName);

long double UniformRand(int lb, int ub);



#endif // __UTILITY_H__