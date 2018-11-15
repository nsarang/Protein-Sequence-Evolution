#ifndef __CPP_LOCATE_H__
#define __CPP_LOCATE_H__


#include <string>

namespace clocate
{


std::string GetExecutablePath();


std::string GetDirectoryPath(const std::string & fullpath);


std::string GetCurrentWorkingDirectory();


bool SetCurrentWorkingDirectory(std::string path);


} // namespace


#endif // __CPP_LOCATE_H__
