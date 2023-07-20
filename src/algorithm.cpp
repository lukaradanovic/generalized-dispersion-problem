#include "algorithm.hpp"
#include <iostream>
#include <algorithm>

Result Algorithm::execute()
{
    startTime = std::chrono::steady_clock::now();
    return multistartVND();
}

int Algorithm::getRunningTime() const
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime).count();
}

bool Algorithm::createGreedySolution(Solution& s)
{
    s.sites.resize(problem.n);
    s.firstCenters.resize(problem.n);
    s.secondCenters.resize(problem.n);
    s.distancesToFirstCenters.resize(problem.n);
    s.distancesToSecondCenters.resize(problem.n);

    std::iota(s.sites.begin(), s.sites.end(), 0);
    
    std::sort(s.sites.begin(), s.sites.end(), 
    [this](int site1, int site2){
        return (problem.capacities[site1] / problem.costs[site1]) > (problem.capacities[site2] / problem.costs[site2]);
    });

    int idx = 0;
    double cost = 0, capacity = 0;
    while (idx < problem.n)
    {
        cost += problem.costs[s.sites[idx]];
        capacity += problem.capacities[s.sites[idx]];
        idx++;
        if (capacity >= problem.minCapacity || cost > problem.maxCost)
            break;
    }
    if (capacity < problem.minCapacity || cost > problem.maxCost)
        return false;

    s.numIncluded = idx;
    s.cost = cost;
    s.capacity = capacity;

    // perform perturbations
    std::uniform_int_distribution<int> distribIn(0, s.numIncluded - 1);
    std::uniform_int_distribution<int> distribOut(s.numIncluded, problem.n - 1);
    for (int i = 0; i < problem.n; i++)
    {
        int idx1 = distribIn(generator), idx2 = distribOut(generator);
        int site1 = s.sites[idx1], site2 = s.sites[idx2];
        if (s.capacity - problem.capacities[site1] + problem.capacities[site2] >= problem.minCapacity && s.cost - problem.costs[site1] + problem.costs[site2] <= problem.maxCost)
        {
            std::swap(s.sites[idx1], s.sites[idx2]);
            s.capacity = s.capacity - problem.capacities[site1] + problem.capacities[site2];
            s.cost = s.cost - problem.costs[site1] + problem.costs[site2];
        }
    }
    
    updateSolutionInfo(problem, s);

    return true;
}

void Algorithm::localSearchSwap(Solution& s, bool& hasImproved)
{
    neighborhoodRuns[0]++;
    int idx1 = -1, idx2 = -1;
    double currentMinDistance = s.minDistance;
    int currentCriticalCount = s.numCritical;
    bool shouldFinish = false;
    for (int i = 0; i < s.numIncluded && !shouldFinish; i++)
    {
        if (!isSiteCritical(s.sites[i], s))
            continue;
        
        for (int j = s.numIncluded; j < problem.n && !shouldFinish; j++)
        {
            if (getRunningTime() > timeMax)
                return;

            if (!isSolutionFeasibleAfterChange({i}, {j}, s, problem))
                continue;

            auto [minDistance, criticalCount] = calcMinDistanceAfterChange({i}, {j}, s, problem);
            
            if (minDistance > currentMinDistance || (minDistance == currentMinDistance && criticalCount < currentCriticalCount))
            {
                currentMinDistance = minDistance;
                currentCriticalCount = criticalCount;
                hasImproved = true;
                idx1 = i;
                idx2 = j;
                shouldFinish = true;
            }
        }
    }

    if (hasImproved)
    {
        std::swap(s.sites[idx1], s.sites[idx2]);
        neighborhoodImprovements[0]++;
        updateSolutionInfo(problem, s);
    }
}

void Algorithm::localSearch2Out1In(Solution& s, bool& hasImproved)
{
    neighborhoodRuns[1]++;
    int idx1 = -1, idx2 = -1, idx3 = -1;
    double currentMinDistance = s.minDistance;
    double currentCriticalCount = s.numCritical;
    bool shouldFinish = false;
    for (int i = 0; i < s.numIncluded && !shouldFinish; i++)
    {
        // at least one site must be critical
        if (!isSiteCritical(s.sites[i], s))
            continue;
        
        for (int j = 0; j < s.numIncluded && !shouldFinish; j++)
        {
            // prevent repeated pairs when both sites are critical 
            if (isSiteCritical(s.sites[j], s) && i >= j)
                continue;

            for (int k = s.numIncluded; k < problem.n && !shouldFinish; k++)
            {
                if (getRunningTime() > timeMax)
                    return;

                if (!isSolutionFeasibleAfterChange({i, j}, {k}, s, problem))
                    continue;

                auto [minDistance, criticalCount] = calcMinDistanceAfterChange({i, j}, {k}, s, problem);
                
                if (minDistance > currentMinDistance || (minDistance == currentMinDistance && criticalCount < currentCriticalCount))
                {
                    currentMinDistance = minDistance;
                    currentCriticalCount = criticalCount;
                    hasImproved = true;
                    idx1 = i;
                    idx2 = j;
                    idx3 = k;
                    shouldFinish = true;
                }
            }
        }
    }

    if (hasImproved)
    {
        std::swap(s.sites[idx1], s.sites[idx3]);
        std::swap(s.sites[idx2], s.sites[s.numIncluded-1]);
        s.numIncluded--;
        neighborhoodImprovements[1]++;
        updateSolutionInfo(problem, s);
    }
}

void Algorithm::localSearch1Out2In(Solution& s, bool& hasImproved)
{
    neighborhoodRuns[2]++;
    int idx1 = -1, idx2 = -1, idx3 = -1;
    double currentMinDistance = s.minDistance;
    double currentCriticalCount = s.numCritical;
    bool shouldFinish = false;
    for (int i = 0; i < s.numIncluded && !shouldFinish; i++)
    {
        if (!isSiteCritical(s.sites[i], s))
            continue;

        for (int j = s.numIncluded; j < problem.n - 1 && !shouldFinish; j++)
        {
            for (int k = j + 1; k < problem.n && !shouldFinish; k++)
            {

                if (getRunningTime() > timeMax)
                    return;

                if (!isSolutionFeasibleAfterChange({i}, {j, k}, s, problem))
                    continue;

                auto [minDistance, criticalCount] = calcMinDistanceAfterChange({i}, {j, k}, s, problem);

                if (minDistance > currentMinDistance || (minDistance == currentMinDistance && criticalCount < currentCriticalCount))
                {
                    currentMinDistance = minDistance;
                    currentCriticalCount = criticalCount;
                    hasImproved = true;
                    idx1 = i;
                    idx2 = j;
                    idx3 = k;
                    shouldFinish = true;
                }
            }
        }
    }

    if (hasImproved)
    {
        std::swap(s.sites[idx1], s.sites[idx2]);
        std::swap(s.sites[s.numIncluded], s.sites[idx3]);
        s.numIncluded++;
        neighborhoodImprovements[2]++;
        updateSolutionInfo(problem, s);
    }
}

void Algorithm::localSearchSwapWithIndexRemembering(Solution& s, bool& hasImproved)
{
    neighborhoodRuns[0]++;
    double currentMinDistance = s.minDistance;
    int currentCriticalCount = s.numCritical;
    int endIndex = s.numIncluded; 
    for (int i = 0; i != endIndex; i++)
    {
        i = i % s.numIncluded;
        if (i == endIndex)
            break;

        if (!isSiteCritical(s.sites[i], s))
            continue;
        
        for (int j = s.numIncluded; j < problem.n; j++)
        {
            if (getRunningTime() > timeMax)
                return;

            if (!isSolutionFeasibleAfterChange({i}, {j}, s, problem))
                continue;

            auto [minDistance, criticalCount] = calcMinDistanceAfterChange({i}, {j}, s, problem);
            
            if (minDistance > currentMinDistance || (minDistance == currentMinDistance && criticalCount < currentCriticalCount))
            {
                currentMinDistance = minDistance;
                currentCriticalCount = criticalCount;

                endIndex = i;
                std::swap(s.sites[i], s.sites[j]);
                neighborhoodImprovements[0]++;
                updateSolutionInfo(problem, s);
                neighborhoodRuns[0]++;
            }
        }
    }
}

void Algorithm::vnd(Solution& s, bool shuffle)
{
    bool hasImproved;
    while (true)
    {
        if (shuffle)
        {
            std::shuffle(s.sites.begin(), std::next(s.sites.begin(), s.numIncluded), generator);
            std::shuffle(std::next(s.sites.begin(), s.numIncluded), s.sites.end(), generator);
        }

        hasImproved = false;
        localSearchSwap(s, hasImproved);
        if (!hasImproved)
            localSearch2Out1In(s, hasImproved);
        if (!hasImproved)
            localSearch1Out2In(s, hasImproved);
        if (!hasImproved)
            break;
        
    }
}

void Algorithm::shake(Solution& s, int k)
{
    for (int i = 0; i < k; i++)
    {
        if (s.numIncluded == problem.n)
            break;

        std::uniform_int_distribution<int> distribOut(s.numIncluded, problem.n - 1);

        int idx = distribOut(generator);
        int site = s.sites[idx];
        if (s.cost + problem.costs[site] <= problem.maxCost)
        {
            std::swap(s.sites[s.numIncluded], s.sites[idx]);
            s.capacity = s.capacity + problem.capacities[site];
            s.cost = s.cost + problem.costs[site];
            s.numIncluded++;
        }
    }
    updateSolutionInfo(problem, s);
}

Result Algorithm::multistartVND()
{
    Result result;

    double bestValue = std::numeric_limits<double>::min();
    Solution bestSol;
    int bestSolTime = -1;
    while (getRunningTime() < timeMax)
    {
        Solution s;
        bool feasibleSolutionCreated = createGreedySolution(s);
        if (!feasibleSolutionCreated)
        {
            if (verbose)
                std::cerr << "Greedy algorithm failed to construct feasible solution." << std::endl;
            return getNonFeasibleSolutionResult();
        }
        vnd(s, true);
        if (s.minDistance > bestValue)
        {
            bestValue = s.minDistance;
            bestSol = s;
            bestSolTime = getRunningTime();

            result.history.push_back({bestValue, bestSolTime});

            if (verbose)
                std::cout << bestValue << "  " << bestSolTime <<  std::endl;
        }
    }

    if (verbose)
    {
        std::cout << "distance: " << bestSol.minDistance << "  cost:  " << bestSol.cost << "  capacity:  " << bestSol.capacity << " numIncluded: " << bestSol.numIncluded << std::endl;
        std::cout << "time: " << bestSolTime << std::endl;
    }

    fillResultInfo(result, bestSol, problem, timeMax, neighborhoodRuns, neighborhoodImprovements);

    return result;
}

Result Algorithm::VNS()
{
    Result result;

    int kMax = problem.n * 0.1;
    double bestValue = std::numeric_limits<double>::min();
    int bestSolTime = -1;

    Solution bestSol;
    bool feasibleSolutionCreated = createGreedySolution(bestSol); 
    if (!feasibleSolutionCreated)
    {
        if (verbose)
            std::cerr << "Greedy algorithm failed to construct feasible solution." << std::endl;
        return getNonFeasibleSolutionResult();
    }
    int k = 1;
    while (k <= kMax && getRunningTime() < timeMax)
    {
        Solution newSol = bestSol;
        shake(newSol, k);
        vnd(newSol);

        if (newSol.minDistance > bestSol.minDistance || (newSol.minDistance == bestSol.minDistance && newSol.numCritical < bestSol.numCritical))
        {
            bestSol = newSol;
            k = 1;

            if (newSol.minDistance > bestValue)
            {
                bestValue = newSol.minDistance;
                bestSol = newSol;
                bestSolTime = getRunningTime();

                result.history.push_back({bestValue, bestSolTime});

                if (verbose)
                    std::cout << bestValue << "  " << bestSolTime <<  std::endl;
            }
        }
        else
        {
            k++;
            if (k > kMax)
                k = 1;
        }
    }

    if (verbose)
    {
        std::cout << "distance: " << bestSol.minDistance << "  cost:  " << bestSol.cost << "  capacity:  " << bestSol.capacity << " numIncluded: " << bestSol.numIncluded << std::endl;
        std::cout << "time: " << bestSolTime << std::endl;
    }

    fillResultInfo(result, bestSol, problem, timeMax, neighborhoodRuns, neighborhoodImprovements);

    return result;
}