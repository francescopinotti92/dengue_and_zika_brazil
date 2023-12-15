//
//  utils.cpp
//  CrossFlavivirus
//
//  Created by user on 31/07/2021.
//

#include "utils.hpp"

// Checks if a given directory exists
bool exists_path (const std::string& name)
{
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

// Creates output (nested) directories and returns path to the innermost directory
std::string path_creator(std::string basic_path, std::vector<std::pair<std::string, int> > folder_info)
{
    std::string mydir = basic_path;
    if (exists_path(mydir)==0) mkdir(mydir.c_str(), ACCESSPERMS);
    for (auto &element :  folder_info)
    {
        mydir += "/" + element.first + "_" + std::to_string(element.second);
        if (exists_path(mydir)==0) mkdir(mydir.c_str(), ACCESSPERMS);
    }
    return(mydir);
}
