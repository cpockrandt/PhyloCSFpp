#pragma once

#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>

void print_info_msg(std::string && msg)
{
    printf("\033[1;34mINFO: %s\033[0m\n", msg.c_str());
//    std::cerr << "\033[1;31m" << "ERROR: " << msg << "\033[0m" << '\n';
}

void print_error_msg(std::string && msg)
{
    printf("\033[1;31mERROR: %s\033[0m\n", msg.c_str());
//    std::cerr << "\033[1;31m" << "ERROR: " << msg << "\033[0m" << '\n';
}

bool create_directory(const std::string & path)
{
    struct stat st = {0};

    if (stat(path.c_str(), &st) == -1)
    {
        mkdir(path.c_str(), 0764);
        return true;
    }
    return false;
}