#include "utility.h"

namespace sp = subprocess;

namespace utility
{

std::vector<std::string> CATH_ListFiles(std::string sDB_Path)
{
    std::vector<std::string> vecProcessedDB, vecBL;
    DIR *hDir;
    dirent *hFile;
    assert(hDir = opendir(sDB_Path.c_str()));
    int total = CountFilesInDir(sDB_Path);
    vecProcessedDB.reserve(total);

    if (FileExists(fl_CATH_BL))
    {
        std::ifstream blFile(fl_CATH_BL);
        std::string sName;
        while (blFile >> sName)
            vecBL.push_back(sName);
    }

    while ((hFile = readdir(hDir)))
    {
        std::string fName = hFile->d_name;
        if (fName.size() != 7 || !isdigit(fName[0]) ||
            std::binary_search(vecBL.begin(), vecBL.end(), fName))
            continue;

        vecProcessedDB.push_back(sDB_Path + fName);
    }
    return vecProcessedDB;
}

void Progress_Indicator(std::string text, long long current, long long total)
{
    std::cout << (current ? "\r\033[F" : "") << text << ": ";
    int percent = current * 100 / total;
    if (percent == 100)
        std::cout << BOLDGREEN << "OK " << RESET << std::endl;
    else
        std::cout << BOLDRED << percent << "%" << RESET << std::endl;
}

/*
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
*/

std::vector<std::string> split(const std::string input, const std::string regex)
{
    // passing -1 as the submatch index parameter performs splitting
    std::regex re(regex);
    std::sregex_token_iterator
        first{input.begin(), input.end(), re, -1},
        last;
    return {first, last};
}

bool startswith(const std::string source, const std::string query)
{
    return source.rfind(query, 0) == 0;
}

bool endswith(const std::string source, const std::string query)
{
    return source.find(query, source.size() - query.size() - 1) != std::string::npos;
}

inline void ltrim(std::string &s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::not1(std::ptr_fun<int, int>(std::isspace))));
}

inline void rtrim(std::string &s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace)))
                .base(),
            s.end());
}

std::string trim(std::string s)
{
    ltrim(s);
    rtrim(s);
    return s;
}

/******************************************************************************
 * Checks to see if a directory exists. Note: This method only checks the
 * existence of the full path AND if path leaf is a dir.
 *
 * @return  >0 if dir exists AND is a dir,
 *           0 if dir does not exist OR exists but not a dir,
 *          <0 if an error occurred like dir isn't accessible (errno is also set)
 *****************************************************************************/
int DirectoryExists(std::string fDir)
{
    struct stat info;

    int statRC = stat(fDir.c_str(), &info);
    if (statRC != 0)
    {
        if (errno == ENOENT)
        {
            return 0;
        } // something along the path does not exist
        if (errno == ENOTDIR)
        {
            return 0;
        } // something in path prefix is not a dir
        return -1;
    }
    return (info.st_mode & S_IFDIR) ? 1 : 0;
}

bool FileExists(std::string fName)
{
    std::ifstream inFile(fName);
    return inFile.good();
}

size_t FileSize(std::istream &isObj)
{
    auto current = isObj.tellg();
    isObj.seekg(0, std::ios::end);
    size_t length = isObj.tellg();
    isObj.seekg(current);
    return length;
}

std::string FileBasename(std::string filename)
{
    size_t last_slash_idx = filename.find_last_of("\\/");
    if (std::string::npos != last_slash_idx)
        filename.erase(0, last_slash_idx + 1);

    size_t period_idx = filename.rfind('.');
    if (std::string::npos != period_idx)
        filename.erase(period_idx);
    return filename;
}

std::string File_md5(std::string fName)
{
#if defined(__linux__)
    std::string ret = sp::check_output({"md5sum", fName.c_str()}).buf.data();
#else
    std::string ret = sp::check_output({"md5", "-r", fName.c_str()}).buf.data();
#endif
    return ret.substr(0, ret.find(' '));
}

int CountFilesInDir(std::string fDir)
{
    auto ls = sp::Popen({"ls", fDir.c_str()}, sp::output{sp::PIPE});
    auto wc = sp::Popen({"wc", "-l"},
                        sp::input{ls.output()},
                        sp::output{sp::PIPE});
    auto res = wc.communicate().first;
    return std::stoi(res.buf.data());
}

long double UniformRand(int lb, int ub)
{
    return lb + (long double)rand() / (RAND_MAX) * (ub - lb);
}

} // namespace utility