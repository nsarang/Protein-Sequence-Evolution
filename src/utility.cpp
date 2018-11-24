#include "utility.h"



std::vector<std::string> CATH_ListFiles(std::string sDB_Path) {
    std::vector<std::string> vecProcessedDB;
    DIR *hDir;
    dirent *hFile;
    assert(hDir = opendir(sDB_Path.c_str()));
    int total = stoi(system_call("ls " + sDB_Path + " | wc -l"));

    vecProcessedDB.reserve(total);
    while ((hFile = readdir(hDir))) {
        std::string fName = hFile->d_name;
        if (fName.size() != 7 || !isdigit(fName[0]))
            continue;

        vecProcessedDB.push_back(sDB_Path + fName);
    }
    return vecProcessedDB;
}


void Progress_Indicator(std::string text, long long current, long long total) {
    if (current == 0)
        std::cout << "|\n";
    std::cout << "\r\033[F" << text << ": ";
    int percent = current * 100 / total;
    if (percent == 100)
        std::cout << BOLDGREEN << "OK " << RESET << std::endl;
    else
        std::cout << BOLDRED << percent << "%" << RESET << std::endl;
}


int system_call_err(std::string command, std::string& stdout) {
    std::array<char, BUFFERSIZE> buffer;

    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        throw std::runtime_error(std::string() + "popen() failed!:\t" + strerror(errno) + "\ncommand: " + command);
    }
    try {
        while (!feof(pipe))
            if (fgets(buffer.data(), BUFFERSIZE, pipe) != NULL)
                stdout += buffer.data();
    } catch (...) {
        pclose(pipe);
        throw;
    }
    int wstat = pclose(pipe);
    pclose(pipe);
    return WEXITSTATUS(wstat);
}


std::string system_call(std::string command) {
    std::array<char, BUFFERSIZE> buffer;
    std::string result = "";

    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
        throw std::runtime_error(std::string() + "popen() failed!:\t" + strerror(errno) + "\ncommand: " + command);
    }
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
#if defined(__linux__)
    std::string ret = subprocess::check_output({"md5sum", fName.c_str()}).buf.data();
#else
    std::string ret = subprocess::check_output({"md5", "-r", fName.c_str()}).buf.data();
#endif
    return ret.substr(0, ret.find(' '));
}


long double UniformRand(int lb, int ub) {
    return lb + (long double)rand() / (RAND_MAX) * (ub - lb);
}