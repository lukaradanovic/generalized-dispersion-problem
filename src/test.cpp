#include "test.hpp"
#include "util.hpp"
#include "structures.hpp"
#include "algorithm.hpp"
#include <random>
#include <filesystem>
#include <iostream>

void processDirectory(const std::string& inputDir, const std::string& outputDir, int testRepetitions)
{
    checkCondition(std::filesystem::create_directory(outputDir), "Couldn't create output directory.");

    std::ofstream outputCSV(outputDir / std::filesystem::path(outputDir + ".csv"));
    checkCondition(!!outputCSV, "Couldn't open CSV file for writing.");
    writeCSVHeader(outputCSV);

    std::random_device rd;
    std::mt19937 generator(rd());

    int counter = 0;
    for (auto & entry : std::filesystem::directory_iterator(inputDir))
    {
        std::cout << ++counter << " " << entry.path().stem() << " ";
        std::filesystem::path inputFilePath = entry.path();

        Problem p = loadData(inputFilePath);

        std::filesystem::path instanceOutputDirPath = outputDir / inputFilePath.stem();
        checkCondition(std::filesystem::create_directory(instanceOutputDirPath), "Couldn't create output directory.");

        for (int i = 1; i <= testRepetitions; i++)
        {
            std::cout << "." << std::flush; // to show progress

            generator.seed(i);

            Algorithm alg = Algorithm(p, generator); 
            Result res = alg.execute();

            std::string outputFilePath = std::string(instanceOutputDirPath / inputFilePath.stem()) + "_result_" + std::to_string(i) + ".txt";
            std::ofstream outputFile(outputFilePath);
            checkCondition(!!outputFile, "Couldn't open file for writing.");

            writeResultToFile(res, outputFile);

            writeResultToCSV(res, outputCSV, std::string(inputFilePath.stem()));
        }

        std::cout << " done." << std::endl;
    }
}

void processFile(const std::string& inputFile, int seed)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    Problem p = loadData(inputFile);
    generator.seed(seed);
    Algorithm alg = Algorithm(p, generator); 
    Result res = alg.execute();
}

void iraceTest(char** argv)
{
    int seed = std::stoi(argv[3]);
    std::string filename(argv[4]);
    double kstepCoef = std::stod(argv[6]);
    double kmaxCoef = std::stod(argv[8]);
    double ls3Prob = std::stod(argv[10]);
    double ls2Prob = std::stod(argv[12]);

    std::random_device rd;
    std::mt19937 generator(rd());
    generator.seed(seed);

    Problem p = loadData(filename);
    Algorithm alg = Algorithm(p, generator, 60000, false, kstepCoef, kmaxCoef, ls3Prob, ls2Prob);
    Result res = alg.execute();

    std::cout << -res.minDistance << std::endl;
}