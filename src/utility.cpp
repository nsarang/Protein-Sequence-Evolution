#include "utility.h"
#include <iostream>
#include <fstream>
#include <array>
#include <algorithm>
#include "cpplocate.h"





double dist(std::tuple<double, double, double> &t1, std::tuple<double, double, double> &t2) {
    double dx = std::get<0>(t1) - std::get<0>(t2),
           dy = std::get<1>(t1) - std::get<1>(t2),
           dz = std::get<2>(t1) - std::get<2>(t2);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}


bool IsStandardAA(std::string abrv)
{
    return !(abrv == "ASX" || abrv == "GLX" || abrv == "SEC" || abrv == "PYL" || abrv == "UNK");
}

/*
void ParsePDB(std::string fPath, std::vector<AminoAcid> &vecAminoAcid) {
    vecAminoAcid.clear();
    std::ifstream inFile(fPath);
    std::string line;

     while (getline(inFile, line)) {
        if (line.size() < 47 || line.substr(0, 4) != "ATOM")
            continue;

        std::string atom_name = trim(line.substr(12, 4));
        if (atom_name != "CA" || !(line[16] == 'A' || line[16] == ' ' || line[16] == '1'))
            continue;

        vecAminoAcid.emplace_back(AminoAcid(
                                         line.substr(17, 3),
                                         std::stod(trim(line.substr(30, 8))),
                                         std::stod(trim(line.substr(38, 8))),
                                         std::stod(trim(line.substr(46, 8)))
                                     ));
        assert(residue_map.count(vecAminoAcid.back().name) > 0);
    }
}
*/

void Progress_Indicator(std::string text, long long current, long long total) {
    std::cout << (current != 0 ?  "\r\033[F" : "" ) << text << ": ";
    int percent = current * 100 / total;
    if (percent == 100)
        std::cout << BOLDGREEN << "OK " << RESET << std::endl;
    else
        std::cout << BOLDRED << percent << "%" << RESET << std::endl;
}


int system_call_err(std::string command, std::string& stdout) {
    std::array<char, BUFFERSIZE> buffer;

    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) throw std::runtime_error("popen() failed!\ncommand: " + command);
    try {
        while (!feof(pipe))
            if (fgets(buffer.data(), BUFFERSIZE, pipe) != NULL)
                stdout += buffer.data();
    } catch (...) {
        pclose(pipe);
        throw;
    }
    int wstat = pclose(pipe);
    return WEXITSTATUS(wstat);
}


std::string system_call(std::string command) {
    std::array<char, BUFFERSIZE> buffer;
    std::string result = "";

    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) throw std::runtime_error("popen() failed!\ncommand: " + command);
    try {
        while (!feof(pipe))
            if (fgets(buffer.data(), BUFFERSIZE, pipe) != NULL)
                result += buffer.data();
    } catch (...) {
        pclose(pipe);
        throw;
    }
    return result;
}


inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
}


inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}


std::string trim(std::string s) {
    ltrim(s);
    rtrim(s);
    return s;
}


bool FileExists(std::string fName) {
    std::ifstream inFile(fName);
    return inFile.good();
}


size_t FileSize(std::istream &isObj) {
    auto current = isObj.tellg();
    isObj.seekg(0, std::ios::end);
    size_t length = isObj.tellg();
    isObj.seekg(current);
    return length;
}


std::string FileBasename(std::string filename) {
    size_t last_slash_idx = filename.find_last_of("\\/");
    if (std::string::npos != last_slash_idx)
        filename.erase(0, last_slash_idx + 1);

    size_t period_idx = filename.rfind('.');
    if (std::string::npos != period_idx)
        filename.erase(period_idx);
    return filename;
}


std::string File_md5(std::string fName) {
    std::string md5 = "md5 -r ";
#if defined(__linux__)
        md5 = "md5sum";
#endif
    std::string ret = system_call("md5 -r " + fName);
    return ret.substr(0, ret.find(' '));
}


long double UniformRand(int lb, int ub) {
    return lb + (long double)rand() / (RAND_MAX) * (ub - lb);
}