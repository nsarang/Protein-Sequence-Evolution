#ifndef __CPP_LOCATE_CPP__
#define __CPP_LOCATE_CPP__


#if defined(_WIN32)
#include <Windows.h>
#elif defined(__linux__)
#include <unistd.h>
#include <limits.h>
#include <dlfcn.h>
#elif defined(__APPLE__)
#include <unistd.h>
#include <limits.h>
#include <mach-o/dyld.h>
#include <dlfcn.h>
#endif

#include <array>
#include <cstdlib>
#include <vector>
#include <iostream>


namespace clocate
{


std::string GetExecutablePath()
{
#if defined(_WIN32)

    std::array<char, MAX_PATH> exePath;

    if (GetModuleFileNameA(GetModuleHandleA(nullptr), exePath.data(), exePath.size()) == 0)
    {
        return "";
    }

    return std::string(exePath.data());


#elif defined(__linux__)

    std::array<char, PATH_MAX> exePath;

    auto len = ::readlink("/proc/self/exe", exePath.data(), exePath.size());

    if (len == -1 || len == exePath.size())
    {
        return "";
    }

    return std::string(exePath.data(), len);


#elif defined(__APPLE__)

    std::array<char, PATH_MAX> exePath;

    auto len = std::uint32_t(exePath.size());

    std::string finalPath;
    if (_NSGetExecutablePath(exePath.data(), &len) == 0)
    {
        auto retPath = realpath(exePath.data(), nullptr);

        if (retPath)
        {
            finalPath = std::string(retPath);
            free(retPath);
        }

    }

    return finalPath;


#else

    return "";

#endif
}


std::string GetDirectoryPath(const std::string & fullpath)
{
    if (fullpath.empty())
    {
        return "";
    }

    auto pos           = fullpath.rfind("/");
    const auto posBack = fullpath.rfind("\\");

    if (pos == std::string::npos || (posBack != std::string::npos && posBack > pos))
    {
        pos = posBack;
    }

    return fullpath.substr(0, pos);
}


std::string GetCurrentWorkingDirectory()
{
#if defined(__linux__) || defined(__APPLE__)

    std::array<char, PATH_MAX> cwdPath;

    if (getcwd(cwdPath.data(), PATH_MAX) != NULL)
    {
        return std::string(cwdPath.data());
    }

    return "";


#elif defined(_WIN32)

    std::array<wchar_t, PATH_MAX> cwdPath;
    if (GetCurrentDirectoryW(MAX_PATH, cwdPath.data()) != 0)
    {
        return std::wstring(cwdPath.data);
    }

    return "";


#else

    return "";

#endif
}


bool SetCurrentWorkingDirectory(std::string path)
{
#if defined(__linux__) || defined(__APPLE__)

    return chdir(path.c_str()) == 0;


#elif defined(_WIN32)

    return SetCurrentDirectory(path.c_str()) != 0;

#else

    return false;

#endif
}


} // namespace


#endif // __CPP_LOCATE_CPP__
