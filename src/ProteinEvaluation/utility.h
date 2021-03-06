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
#include <mutex>
#include <sys/types.h>
#include <sys/stat.h>
#include <regex>
#include "constants.h"
#include "cpplocate.hpp"
#include "mapping.h"
#include "subprocess.hpp"


namespace utility {

const char RESET[] = "\033[0m";
const char BOLDRED[] = "\033[1m\033[31m";
const char BOLDGREEN[] = "\033[1m\033[32m";
const int BUFFERSIZE = 1024;


template<class vecType>
std::vector<vecType> VectorSubset(const std::vector<vecType>& vec,
    							  const std::vector<size_t>& indices);


template<class FuncType, class vecType>
void Thread_Manager(std::vector<std::function<FuncType> > vecFunctions,
	                std::vector<vecType>& vecDatabase,
	                bool bVerbose = false, std::string sMsg = "",
	                int nThreads = std::thread::hardware_concurrency());


template<class FuncType, class vecType>
void Processing_Thread(std::vector<std::function<FuncType> > vecFunctions,
	                   std::vector<vecType> vecBatch, int& nCount_Now);


std::vector<std::string> CATH_ListFiles(std::string path);

void Progress_Indicator(std::string text, long long current, long long total);

// int system_call_err(std::string command, std::string& stdout);

std::vector<std::string> split(const std::string input, const std::string regex);

bool startswith(const std::string source, const std::string query);

bool endswith(const std::string source, const std::string query);

void ltrim(std::string &s);

void rtrim(std::string &s);

std::string trim(std::string s);

int DirectoryExists(std::string fDir);

bool FileExists(std::string fName);

size_t FileSize(std::istream &isObj);

std::string FileBasename(std::string filename);

std::string File_md5(std::string fName);

int CountFilesInDir(std::string fDir);

long double UniformRand(int lb, int ub);


} // namespace


#include "utility.tpp"

#endif // __UTILITY_H__