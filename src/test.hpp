#ifndef _TEST_H_
#define _TEST_H_ 

#include <string>

void processDirectory(const std::string& inputDir, const std::string& outputDir, int testRepetitions);
void processFile(const std::string& inputFile, int seed);
void iraceTest(char** argv);
void processDirectoryWithParams(const std::string& inputDir, const std::string& outputDir, int testRepetitions);

#endif