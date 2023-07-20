#include "util.hpp"
#include <iostream>

void checkCondition(bool cond, const std::string& failMessage)
{
    if (!cond)
    {
        std::cerr << failMessage << std::endl;
        exit(EXIT_FAILURE);
    }
}