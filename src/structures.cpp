#include "structures.hpp"
#include "util.hpp"
#include <limits>
#include <algorithm>

Problem loadData(const std::string& filename)
{
    Problem problem;

    std::ifstream input(filename);
    checkCondition(!!input, "Opening file failed.");

    int n;
    input >> n;
    problem.n = n;
    problem.costs.resize(n);
    problem.capacities.resize(n);
    problem.distances.resize(n);
    for (int i = 0; i < n; i++)
    {
        problem.distances[i].resize(n);
        problem.distances[i][i] = std::numeric_limits<double>::max();
    }

    int numPairs = (n * (n - 1)) / 2, idx1, idx2;
    double distance;
    for (int i = 0; i < numPairs; i++)
    {
        input >> idx1 >> idx2 >> distance;

        // indices in input format start from 1
        problem.distances[idx1 - 1][idx2 - 1] = distance;
        problem.distances[idx2 - 1][idx1 - 1] = distance; 
    }

    double placeholder; // for reading values not needed for this problem
    for (int i = 0; i < n; i++)
    {
        input >> placeholder >> problem.costs[i] >> placeholder >> problem.capacities[i];
    }
    input >> problem.maxCost >> placeholder >> problem.minCapacity;

    return problem;
}

void updateSolutionInfo(const Problem& p, Solution& s)
{
    for (int i = 0; i < p.n; i++)
    {
        double distToFirst = std::numeric_limits<double>::max(), distToSecond = std::numeric_limits<double>::max();
        int firstCenter = -1, secondCenter = -1; 

        for (int j = 0; j < s.numIncluded; j++)
        {
            double dist = p.distances[s.sites[i]][s.sites[j]];
            if (dist < distToFirst)
            {
                distToSecond = distToFirst;
                secondCenter = firstCenter;
                
                distToFirst = dist;
                firstCenter = s.sites[j];
            }
            else if (dist < distToSecond)
            {
                distToSecond = dist;
                secondCenter = s.sites[j];
            }
        }

        s.firstCenters[s.sites[i]] = firstCenter;
        s.secondCenters[s.sites[i]] = secondCenter;
        s.distancesToFirstCenters[s.sites[i]] = distToFirst;
        s.distancesToSecondCenters[s.sites[i]] = distToSecond;
    }

    s.minDistance = std::numeric_limits<double>::max();
    s.cost = 0;
    s.capacity = 0;
    for (int i = 0; i < s.numIncluded; i++)
    {
        s.minDistance = std::min(s.minDistance, s.distancesToFirstCenters[s.sites[i]]);
        s.cost += p.costs[s.sites[i]];
        s.capacity += p.capacities[s.sites[i]];
    }

    s.numCritical = 0;
    for (int i = 0; i < s.numIncluded; i++)
        s.numCritical += s.distancesToFirstCenters[s.sites[i]] == s.minDistance ? 1 : 0;
}

bool isSiteCritical(int site, const Solution& s)
{
    return s.distancesToFirstCenters[site] == s.minDistance;
}

DistanceAndCriticalCount calcMinDistanceAfterChange(const std::vector<int>& removedIndices, const std::vector<int>& addedIndices, const Solution& s, const Problem& p)
{
    double minDistance = std::numeric_limits<double>::max();
    int criticalCount = 0;

    std::vector<int> removedSites(removedIndices);
    std::transform(removedSites.begin(), removedSites.end(), removedSites.begin(), [&s](int idx){ return s.sites[idx];});
    for (int idx = 0; idx < s.numIncluded; idx++)
    {
        if (std::find(removedIndices.begin(), removedIndices.end(), idx) != removedIndices.end())
            continue;
        
        double dist = calcChangedDistanceForOneSite(s.sites[idx], removedIndices, removedSites, addedIndices, s, p);
        updateTempMinDistanceAndCriticalCount(dist, minDistance, criticalCount);
    }

    for (int idx : addedIndices)
    {
        double dist = calcChangedDistanceForOneSite(s.sites[idx], removedIndices, removedSites, addedIndices, s, p);
        updateTempMinDistanceAndCriticalCount(dist, minDistance, criticalCount);
    }

    return {minDistance, criticalCount};
}

void updateTempMinDistanceAndCriticalCount(double newDist, double& minDistance, int& criticalCount)
{
    if (newDist < minDistance)
    {
        minDistance = newDist;
        criticalCount = 1;
    }
    else if (newDist == minDistance)
    {
        criticalCount++;
    }
}

double calcChangedDistanceForOneSite(int site, const std::vector<int>& removedIndices, const std::vector<int>& removedSites, const std::vector<int>& addedIndices, const Solution& s, const Problem& p)
{
    double dist = std::numeric_limits<double>::max();
    if (std::find(removedSites.begin(), removedSites.end(), s.firstCenters[site]) == removedSites.end())
    {
        dist = s.distancesToFirstCenters[site];
        for (int addedIdx : addedIndices)
            dist = std::min(dist, p.distances[site][s.sites[addedIdx]]);
    }
    else if (std::find(removedSites.begin(), removedSites.end(), s.secondCenters[site]) == removedSites.end()) {
        dist = s.distancesToSecondCenters[site];
        for (int addedIdx : addedIndices)
            dist = std::min(dist, p.distances[site][s.sites[addedIdx]]);
    }
    else
    {
        for (int l = 0; l < s.numIncluded; l++)
        {
            if (std::find(removedIndices.begin(), removedIndices.end(), l) != removedIndices.end())
                continue;
            dist = std::min(dist, p.distances[site][s.sites[l]]);
        }
        for (int addedIdx : addedIndices)
            dist = std::min(dist, p.distances[site][s.sites[addedIdx]]);
    }
    return dist;
}

bool isSolutionFeasibleAfterChange(const std::vector<int>& removedIndices, const std::vector<int>& addedIndices, const Solution& s, const Problem& p)
{
    double capacity = s.capacity, cost = s.cost;
    for (int idx : removedIndices)
    {
        capacity -= p.capacities[s.sites[idx]];
        cost -= p.costs[s.sites[idx]];
    }
    for (int idx : addedIndices)
    {
        capacity += p.capacities[s.sites[idx]];
        cost += p.costs[s.sites[idx]];
    }
    return capacity >= p.minCapacity && cost <= p.maxCost;
}

std::vector<int> getIncludedSitesForOutput(const Solution& s)
{
    std::vector<int> included;
    for (int i = 0; i < s.numIncluded; i++)
        included.push_back(s.sites[i] + 1); // labeling starts from 1
    return included;
}

void fillResultInfo(Result& result, const Solution& bestSol, const Problem& p, int timeMax, const std::vector<int>& nRuns, const std::vector<int>& nImps)
{
    result.cost = bestSol.cost;
    result.capacity = bestSol.capacity;
    result.time = result.history.size() > 0 ? result.history[result.history.size() - 1].second : -1;
    result.minDistance = bestSol.minDistance;
    result.maxTime = timeMax;
    result.minCapacity = p.minCapacity;
    result.maxCost = p.maxCost;
    result.includedSites = getIncludedSitesForOutput(bestSol);
    result.neighborhoodImprovements = nImps;
    result.neighborhoodRuns = nRuns;
}

Result getNonFeasibleSolutionResult()
{
    Result result;
    result.minDistance = -1;
    result.cost = -1;
    result.capacity = -1;
    result.time = -1;
    result.minDistance = -1;
    result.maxTime = -1;
    result.minCapacity = -1;
    result.maxCost = -1;
    result.neighborhoodImprovements = {0, 0, 0};
    result.neighborhoodRuns = {0, 0, 0};
    return result;
}

void writeResultToCSV(const Result& res, std::ofstream& output, const std::string& instanceName)
{
    output << instanceName << "," << res.minDistance << "," << res.time << "," << res.maxTime << "," << res.cost << "," << res.maxCost << "," << res.capacity << "," << res.minCapacity 
    << "," << res.neighborhoodImprovements[0] <<"/"<< res.neighborhoodRuns[0]<<","<< res.neighborhoodImprovements[1] <<"/"<< res.neighborhoodRuns[1]<<","<< res.neighborhoodImprovements[2] <<"/"<< res.neighborhoodRuns[2]<< std::endl;
}

void writeCSVHeader(std::ofstream& output)
{
    output << "Instance, Min Distance, Time To Best (ms), Max Time (ms), Cost, Max Cost, Capacity, Min Capacity, LS1, LS2, LS3" << std::endl;
}

void writeResultToFile(const Result& result, std::ofstream& output)
{   
    for (const DistanceAndTime & entry : result.history)
        output << entry.first << ":" << entry.second << " ";
    output << std::endl;

    for (int site : result.includedSites)
        output << site << " ";
    output << std::endl;

    output << "minDistance: " << result.minDistance << std::endl;
    output << "cost: " << result.cost << std::endl;
    output << "capacity: " << result.capacity << std::endl;
    output << "time: " << result.time << std::endl;
    output << "maxTime: " << result.maxTime << std::endl;
    output << "LS1: " << result.neighborhoodImprovements[0] << "/" << result.neighborhoodRuns[0] << std::endl;
    output << "LS2: " << result.neighborhoodImprovements[1] << "/" << result.neighborhoodRuns[1] << std::endl;
    output << "LS3: " << result.neighborhoodImprovements[2] << "/" << result.neighborhoodRuns[2] << std::endl;
}