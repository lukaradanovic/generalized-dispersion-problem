#include "util.hpp"
#include "test.hpp"
#include <string>

int main(int argc, char **argv) 
{
    checkCondition(argc >= 2, "Input file path missing.");
    std::string filename(argv[1]);

    processDirectory(filename, "multistartVND_noLS3", 3);
    //processFile(filename, 4);

    return 0;
}