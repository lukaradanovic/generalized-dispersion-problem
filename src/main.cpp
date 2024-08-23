#include "util.hpp"
#include "test.hpp"
#include <string>

int main(int argc, char **argv) 
{
    checkCondition(argc >= 2, "Input file path missing.");
    std::string filename(argv[1]);

    processDirectory(filename, "/out/VNS", 30);
	//processDirectoryWithParams(filename, "/out/VNS", 10);
    //processFile(filename, 1);

    /*Problem p = loadData(argv[1]);
    std::cout << getMinDistance(readIncludedSites(argv[2]), p) << std::endl;*/

    return 0;
}