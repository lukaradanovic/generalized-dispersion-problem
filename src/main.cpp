#include "util.hpp"
#include "test.hpp"
#include <string>

int main(int argc, char **argv) 
{
    checkCondition(argc >= 2, "Input file path missing.");
    std::string filename(argv[1]);

    processDirectory(filename, "VNS_kstep025_kmax01_shake05_075_noLS3", 3);
    //processFile(filename, 4);

    return 0;
}