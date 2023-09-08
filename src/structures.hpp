#ifndef _STRUCTURES_H_
#define _STRUCTURES_H_ 

#include <utility>
#include <fstream>
#include <vector>
#include <string>

using DistanceAndCriticalCount = std::pair<double, int>;
using DistanceAndTime = std::pair<double, int>;

struct Problem
{
    std::vector<std::vector<double>> distances;
    std::vector<double> capacities;
    std::vector<double> costs;
    double minCapacity;
    double maxCost;
    int n;
};

struct Solution
{
    std::vector<int> sites;
    std::vector<int> firstCenters;
    std::vector<int> secondCenters;
    std::vector<double> distancesToFirstCenters;
    std::vector<double> distancesToSecondCenters;
    int numIncluded;
    int numCritical;
    double cost;
    double capacity;
    double minDistance;
};

struct Result
{
    std::vector<int> includedSites;
    std::vector<DistanceAndTime> history;
    std::vector<int> neighborhoodRuns;
    std::vector<int> neighborhoodImprovements;
    double minDistance;
    double cost;
    double capacity;
    double time;
    double maxCost;
    double minCapacity;
    double maxTime;
    int shakeIters;
    int shakeLoops;
};

Problem loadData(const std::string& filename);

void updateSolutionInfo(const Problem& p, Solution& s);
bool isSiteCritical(int site, const Solution& s);
void updateTempMinDistanceAndCriticalCount(double newDist, double& minDistance, int& criticalCount);
bool isSolutionFeasibleAfterChange(const std::vector<int>& removedIndices, const std::vector<int>& addedIndices, const Solution& s, const Problem& p);

std::vector<int> getIncludedSitesForOutput(const Solution& s);

void fillResultInfo(Result& result, const Solution& bestSol, const Problem& p, int timeMax, const std::vector<int>& nRuns, const std::vector<int>& nImps, int shakeIters = 0, int shakeLoops = 0);
Result getNonFeasibleSolutionResult();
void writeResultToCSV(const Result& res, std::ofstream& output, const std::string& instanceName);
void writeCSVHeader(std::ofstream& output);
void writeResultToFile(const Result& result, std::ofstream& output);

double getMinDistance(const std::vector<int> includedSites, const Problem& p);
std::vector<int> readIncludedSites(const std::string& filename);

#endif