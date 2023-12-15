//
//  utils.hpp
//  CrossFlavivirus
//
//  Created by MacBook Pro on 06/01/2021.
//

#ifndef utils_h
#define utils_h

#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>


// Checks if a given directory exists
bool exists_path (const std::string& name);

// Creates output (nested) directories and returns path to the innermost directory
std::string path_creator(std::string basic_path, std::vector<std::pair<std::string, int> > folder_info);



#endif /* utils_h */
