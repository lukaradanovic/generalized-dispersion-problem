#include "util.hpp"
#include "test.hpp"
#include <string>

int main(int argc, char **argv) 
{
    checkCondition(argc >= 2, "Input file path missing.");
    std::string filename(argv[1]);

    processDirectory(filename, "VNS_kstep05_kmax01_shake052_024_noLS3", 3);
    //processFile(filename, 1);

    /*Problem p = loadData(argv[1]);
    std::cout << getMinDistance(readIncludedSites(argv[2]), p) << std::endl;*/

    return 0;
}