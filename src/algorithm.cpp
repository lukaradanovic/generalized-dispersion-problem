#include "algorithm.hpp"
#include <iostream>
#include <algorithm>

Result Algorithm::execute()
{
    startTime = std::chrono::steady_clock::now();
    return VNS();
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
    int outIdx = -1, inIdx = -1;
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

            auto [minDistance, criticalCount] = getSwapEvaluation(s.sites[i], s.sites[j], s);
            
            if (minDistance > currentMinDistance || (minDistance == currentMinDistance && criticalCount < currentCriticalCount))
            {
                currentMinDistance = minDistance;
                currentCriticalCount = criticalCount;
                hasImproved = true;
                outIdx = i;
                inIdx = j;
                shouldFinish = true;
            }
        }
    }

    if (hasImproved)
    {
        performSwap(outIdx, inIdx, s);
        neighborhoodImprovements[0]++;
    }
}

void Algorithm::localSearch2Out1In(Solution& s, bool& hasImproved)
{
    neighborhoodRuns[1]++;
    int out1Idx = -1, out2Idx = -1, inIdx = -1;
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

                auto [minDistance, criticalCount] = get2Out1InEvaluation(s.sites[i], s.sites[j], s.sites[k], s); //calcMinDistanceAfterChange({i, j}, {k}, s, problem);
                
                if (minDistance > currentMinDistance || (minDistance == currentMinDistance && criticalCount < currentCriticalCount))
                {
                    currentMinDistance = minDistance;
                    currentCriticalCount = criticalCount;
                    hasImproved = true;
                    out1Idx = i;
                    out2Idx = j;
                    inIdx = k;
                    shouldFinish = true;
                }
            }
        }
    }

    if (hasImproved)
    {
        perform2Out1In(out1Idx, out2Idx, inIdx, s);
        neighborhoodImprovements[1]++;
    }
}

void Algorithm::localSearch1Out2In(Solution& s, bool& hasImproved)
{
    neighborhoodRuns[2]++;
    int outIdx = -1, in1Idx = -1, in2Idx = -1;
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

                auto [minDistance, criticalCount] = get1Out2InEvaluation(s.sites[i], s.sites[j], s.sites[k], s); //calcMinDistanceAfterChange({i}, {j, k}, s, problem);

                if (minDistance > currentMinDistance || (minDistance == currentMinDistance && criticalCount < currentCriticalCount))
                {
                    currentMinDistance = minDistance;
                    currentCriticalCount = criticalCount;
                    hasImproved = true;
                    outIdx = i;
                    in1Idx = j;
                    in2Idx = k;
                    shouldFinish = true;
                }
            }
        }
    }

    if (hasImproved)
    {
        perform1Out2In(outIdx, in1Idx, in2Idx, s);
        neighborhoodImprovements[2]++;
    }
}

void Algorithm::vnd(Solution& s, bool shuffle, bool doLS3)
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
        if (!hasImproved && doLS3)
            localSearch1Out2In(s, hasImproved);
        if (!hasImproved)
            break;
        
    }
}

void Algorithm::shake(Solution& s, int k)
{
    for (int i = 0; i < k; i++)
    {
        std::shuffle(s.sites.begin(), std::next(s.sites.begin(), s.numIncluded), generator);
        std::shuffle(std::next(s.sites.begin(), s.numIncluded), s.sites.end(), generator);

        std::uniform_real_distribution<> distrib(0.0, 1.0);
        double randNum = distrib(generator);

        bool found = false;

        if (randNum < ls3Prob)
            shake1Out2In(s);
        else if (randNum < ls3Prob + ls2Prob)
            shake2Out1In(s);
        else
            shakeSwap(s);
    }
    // shaking in a neighborhood updates cost and capacity, but information on centers still needs to be updated
    updateSolutionInfo(problem, s);
}

void Algorithm::shakeSwap(Solution& s)
{
    for (int idx1 = 0; idx1 < s.numIncluded; idx1++)
    {
        for (int idx2 = s.numIncluded; idx2 < problem.n; idx2++)
        {
            int site1 = s.sites[idx1], site2 = s.sites[idx2];
            if (s.capacity - problem.capacities[site1] + problem.capacities[site2] >= problem.minCapacity
                && s.cost - problem.costs[site1] + problem.costs[site2] <= problem.maxCost)
            {
                s.capacity = s.capacity - problem.capacities[site1] + problem.capacities[site2];
                s.cost = s.cost - problem.costs[site1] + problem.costs[site2];
                std::swap(s.sites[idx1], s.sites[idx2]);
                return;
            }
        }
    }
}

void Algorithm::shake2Out1In(Solution& s)
{
    for (int idx1 = 0; idx1 < s.numIncluded - 1; idx1++)
    {
        for (int idx2 = idx1 + 1; idx2 < s.numIncluded; idx2++)
        {
            for (int idx3 = s.numIncluded; idx3 < problem.n; idx3++)
            {
                int site1 = s.sites[idx1], site2 = s.sites[idx2], site3 = s.sites[idx3];
                if (s.capacity - problem.capacities[site1] - problem.capacities[site2] + problem.capacities[site3] >= problem.minCapacity
                    && s.cost - problem.costs[site1] - problem.costs[site2] + problem.costs[site3] <= problem.maxCost)
                {
                    s.capacity = s.capacity - problem.capacities[site1] - problem.capacities[site2] + problem.capacities[site3];
                    s.cost = s.cost - problem.costs[site1] - problem.costs[site2] + problem.costs[site3];
                    std::swap(s.sites[idx1], s.sites[idx3]);
                    std::swap(s.sites[idx2], s.sites[s.numIncluded-1]);
                    s.numIncluded--;
                    return;
                }
            }
        }
    }
}

void Algorithm::shake1Out2In(Solution& s)
{
    for (int idx1 = 0; idx1 < s.numIncluded; idx1++)
    {
        for (int idx2 = s.numIncluded; idx2 < problem.n - 1; idx2++)
        {
            for (int idx3 = idx2 + 1; idx3 < problem.n; idx3++)
            {
                int site1 = s.sites[idx1], site2 = s.sites[idx2], site3 = s.sites[idx3];
                if (s.capacity - problem.capacities[site1] + problem.capacities[site2] + problem.capacities[site3] >= problem.minCapacity
                    && s.cost - problem.costs[site1] + problem.costs[site2] + problem.costs[site3] <= problem.maxCost)
                {
                    s.capacity = s.capacity - problem.capacities[site1] + problem.capacities[site2] + problem.capacities[site3];
                    s.cost = s.cost - problem.costs[site1] + problem.costs[site2] + problem.costs[site3];
                    std::swap(s.sites[idx1], s.sites[idx2]);
                    std::swap(s.sites[s.numIncluded], s.sites[idx3]);
                    s.numIncluded++;
                    return;
                }
            }
        }
    }
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
        //addAditionalSites(s, 0.2);
        vnd(s, true, false);
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

    int kMax = 10;//problem.n * kmaxCoef;
    int kMin = 1;
    int kStep = 1;//kMax * kstepCoef;
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
    int shakeIters = 0;
    int shakeLoops = 0;
    int k = kMin;
    while (getRunningTime() < timeMax)
    {
        Solution newSol = bestSol;
        shake(newSol, k);
        vnd(newSol, true, doLS3);

        if (newSol.minDistance > bestSol.minDistance || (newSol.minDistance == bestSol.minDistance && newSol.numCritical < bestSol.numCritical))
        {
            bestSol = newSol;
            k = kMin;

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
            shakeIters++;
            k += kStep;
            if (k > kMax)
            {
                shakeLoops++;
                if (kMin == kStep)
                    kMin = 1;
                else
                    kMin++;
                k = kMin;
            }
        }
    }

    if (verbose)
    {
        std::cout << "distance: " << bestSol.minDistance << "  cost:  " << bestSol.cost << "  capacity:  " << bestSol.capacity << " numIncluded: " << bestSol.numIncluded << std::endl;
        std::cout << "time: " << bestSolTime << std::endl;
    }

    fillResultInfo(result, bestSol, problem, timeMax, neighborhoodRuns, neighborhoodImprovements, shakeIters, shakeLoops);

    return result;
}

DistanceAndCriticalCount Algorithm::getSwapEvaluation(int out, int in, const Solution& s)
{
    double minDist = std::numeric_limits<double>::max();
    int criticalCount = 0;
    double dist;
    for (int i = 0; i < s.numIncluded; i++)
    {
        int site = s.sites[i];
        if (site == out)
            continue;
        if (s.firstCenters[site] != out)
            dist = std::min(s.distancesToFirstCenters[site], problem.distances[site][in]);
        else
            dist = std::min(s.distancesToSecondCenters[site], problem.distances[site][in]);
        
        updateTempMinDistanceAndCriticalCount(dist, minDist, criticalCount);
    }
    if (s.firstCenters[in] != out)
        dist = s.distancesToFirstCenters[in];
    else
        dist = s.distancesToSecondCenters[in];
    updateTempMinDistanceAndCriticalCount(dist, minDist, criticalCount);

    return {minDist, criticalCount};
}

DistanceAndCriticalCount Algorithm::get2Out1InEvaluation(int out1, int out2, int in, const Solution& s)
{
    double minDist = std::numeric_limits<double>::max();
    int criticalCount = 0;
    double dist;
    for (int i = 0; i < s.numIncluded; i++)
    {
        int site = s.sites[i];
        if (site == out1 || site == out2)
            continue;
        if (s.firstCenters[site] != out1 && s.firstCenters[site] != out2)
            dist = std::min(s.distancesToFirstCenters[site], problem.distances[site][in]);
        else if (s.secondCenters[site] != out1 && s.secondCenters[site] != out2)
            dist = std::min(s.distancesToSecondCenters[site], problem.distances[site][in]);
        else
        {
            dist = problem.distances[site][in];
            for (int j = 0; j < s.numIncluded; j++)
            {
                int site2 = s.sites[j];
                if (site2 == out1 || site2 == out2)
                    continue;
                dist = std::min(dist, problem.distances[site][site2]); 
            }
        }

        updateTempMinDistanceAndCriticalCount(dist, minDist, criticalCount);
    }

    if (s.firstCenters[in] != out1 && s.firstCenters[in] != out2)
        dist = s.distancesToFirstCenters[in];
    else if (s.secondCenters[in] != out1 && s.secondCenters[in] != out2)
        dist = s.distancesToSecondCenters[in];
    else
    {
        dist = std::numeric_limits<double>::max();
        for (int j = 0; j < s.numIncluded; j++)
        {
            int site2 = s.sites[j];
            if (site2 == out1 || site2 == out2)
                continue;
            dist = std::min(dist, problem.distances[in][site2]); 
        }
    }

    updateTempMinDistanceAndCriticalCount(dist, minDist, criticalCount);

    return {minDist, criticalCount};
}

DistanceAndCriticalCount Algorithm::get1Out2InEvaluation(int out, int in1, int in2, const Solution& s)
{
    double minDist = std::numeric_limits<double>::max();
    int criticalCount = 0;
    double dist;
    for (int i = 0; i < s.numIncluded; i++)
    {
        int site = s.sites[i];
        if (site == out)
            continue;

        if (s.firstCenters[site] != out)
            dist = std::min(s.distancesToFirstCenters[site], std::min(problem.distances[site][in1], problem.distances[site][in2]));
        else
            dist = std::min(s.distancesToSecondCenters[site], std::min(problem.distances[site][in1], problem.distances[site][in2]));

        updateTempMinDistanceAndCriticalCount(dist, minDist, criticalCount);
    }

    if (s.distancesToFirstCenters[in1] != out)
        dist = std::min(s.distancesToFirstCenters[in1], problem.distances[in1][in2]);
    else
        dist = std::min(s.distancesToSecondCenters[in1], problem.distances[in1][in2]);
    updateTempMinDistanceAndCriticalCount(dist, minDist, criticalCount);

    if (s.distancesToFirstCenters[in2] != out)
        dist = std::min(s.distancesToFirstCenters[in2], problem.distances[in2][in1]);
    else
        dist = std::min(s.distancesToSecondCenters[in2], problem.distances[in2][in1]);
    updateTempMinDistanceAndCriticalCount(dist, minDist, criticalCount);

    return {minDist, criticalCount};
}

void Algorithm::performSwap(int outIdx, int inIdx, Solution& s)
{
    int out = s.sites[outIdx], in = s.sites[inIdx];
    std::swap(s.sites[outIdx], s.sites[inIdx]);
    s.cost = s.cost - problem.costs[out] + problem.costs[in];
    s.capacity = s.capacity - problem.capacities[out] + problem.capacities[in];

    s.minDistance = std::numeric_limits<double>::max();
    s.numCritical = 0;

    for (int i = 0; i < problem.n; i++)
    {
        int site = s.sites[i];
        if (s.firstCenters[site] == out)
        {
            if (problem.distances[site][in] <= s.distancesToSecondCenters[site])
            {
                s.firstCenters[site] = in;
                s.distancesToFirstCenters[site] = problem.distances[site][in]; 
            }
            else
            {
                s.firstCenters[site] = s.secondCenters[site];
                s.distancesToFirstCenters[site] = s.distancesToSecondCenters[site];

                calculateSecondCenter(site, s);
            }
        }
        else
        {
            if (problem.distances[site][in] < s.distancesToFirstCenters[site])
            {
                s.secondCenters[site] = s.firstCenters[site];
                s.distancesToSecondCenters[site] = s.distancesToFirstCenters[site];

                s.firstCenters[site] = in;
                s.distancesToFirstCenters[site] = problem.distances[site][in]; 
            }
            else if (problem.distances[site][in] < s.distancesToSecondCenters[site])
            {
                s.secondCenters[site] = in;
                s.distancesToSecondCenters[site] = problem.distances[site][in]; 
            }
            else if (s.secondCenters[site] == out)
            {
                calculateSecondCenter(site, s);
            }
        }

        if (i < s.numIncluded)
            updateTempMinDistanceAndCriticalCount(s.distancesToFirstCenters[site], s.minDistance, s.numCritical);
    }
}

void Algorithm::perform2Out1In(int out1Idx, int out2Idx, int inIdx, Solution& s)
{
    int out1 = s.sites[out1Idx], out2 = s.sites[out2Idx], in = s.sites[inIdx]; 
    std::swap(s.sites[out1Idx], s.sites[inIdx]);
    std::swap(s.sites[out2Idx], s.sites[s.numIncluded-1]);
    s.numIncluded--;
    s.cost = s.cost - problem.costs[out1] - problem.costs[out2] + problem.costs[in];
    s.capacity = s.capacity - problem.capacities[out1] - problem.capacities[out2] + problem.capacities[in];

    s.minDistance = std::numeric_limits<double>::max();
    s.numCritical = 0;

    for (int i = 0; i < problem.n; i++)
    {
        int site = s.sites[i];
        if (s.firstCenters[site] != out1 && s.firstCenters[site] != out2)
        {
            if (problem.distances[site][in] < s.distancesToFirstCenters[site])
            {
                s.secondCenters[site] = s.firstCenters[site];
                s.distancesToSecondCenters[site] = s.distancesToFirstCenters[site];

                s.firstCenters[site] = in;
                s.distancesToFirstCenters[site] = problem.distances[site][in]; 
            }
            else if (problem.distances[site][in] < s.distancesToSecondCenters[site])
            {
                s.secondCenters[site] = in;
                s.distancesToSecondCenters[site] = problem.distances[site][in]; 
            }
            else if (s.secondCenters[site] == out1 || s.secondCenters[site] == out2)
            {
                calculateSecondCenter(site, s);
            }
        }
        else
        {
            if (s.secondCenters[site] != out1 && s.secondCenters[site] != out2)
            {
                if (problem.distances[site][in] < s.distancesToSecondCenters[site])
                {
                    s.firstCenters[site] = in;
                    s.distancesToFirstCenters[site] = problem.distances[site][in];
                }
                else
                {
                    s.firstCenters[site] = s.secondCenters[site];
                    s.distancesToFirstCenters[site] = s.distancesToSecondCenters[site];

                    calculateSecondCenter(site, s);
                }
            }
            else
            {
                s.distancesToFirstCenters[site] = std::numeric_limits<double>::max();
                s.distancesToSecondCenters[site] = std::numeric_limits<double>::max();
                for (int j = 0; j < s.numIncluded; j++)
                {
                    int site2 = s.sites[j];
                    double dist = problem.distances[site][site2];
                    if (dist < s.distancesToFirstCenters[site])
                    {
                        s.secondCenters[site] = s.firstCenters[site];
                        s.distancesToSecondCenters[site] = s.distancesToFirstCenters[site];

                        s.firstCenters[site] = site2;
                        s.distancesToFirstCenters[site] = dist; 
                    }
                    else if (dist < s.distancesToSecondCenters[site])
                    {
                        s.secondCenters[site] = site2;
                        s.distancesToSecondCenters[site] = dist; 
                    }
                }
            }
        }

        if (i < s.numIncluded)
            updateTempMinDistanceAndCriticalCount(s.distancesToFirstCenters[site], s.minDistance, s.numCritical);
    }
}

void Algorithm::perform1Out2In(int outIdx, int in1Idx, int in2Idx, Solution& s)
{
    // 1Out2In can be performed as a composition of swap and insertion
    performSwap(outIdx, in1Idx, s);

    int in = s.sites[in2Idx];
    std::swap(s.sites[s.numIncluded], s.sites[in2Idx]);
    s.numIncluded++;
    s.cost += problem.costs[in];
    s.capacity += problem.capacities[in];

    s.minDistance = std::numeric_limits<double>::max();
    s.numCritical = 0;

    for (int i = 0; i < problem.n; i++)
    {
        int site = s.sites[i];
        if (problem.distances[site][in] < s.distancesToFirstCenters[site])
        {
            s.secondCenters[site] = s.firstCenters[site];
            s.distancesToSecondCenters[site] = s.distancesToFirstCenters[site];

            s.firstCenters[site] = in;
            s.distancesToFirstCenters[site] = problem.distances[site][in]; 
        }
        else if (problem.distances[site][in] < s.distancesToSecondCenters[site])
        {
            s.secondCenters[site] = in;
            s.distancesToSecondCenters[site] = problem.distances[site][in];
        }

        if (i < s.numIncluded)
            updateTempMinDistanceAndCriticalCount(s.distancesToFirstCenters[site], s.minDistance, s.numCritical);
    }
}

void Algorithm::calculateSecondCenter(int site, Solution& s)
{
    s.distancesToSecondCenters[site] = std::numeric_limits<double>::max();
    for (int j = 0; j < s.numIncluded; j++)
    {
        int site2 = s.sites[j];
        if (site2 == s.firstCenters[site])
            continue;
        if (problem.distances[site][site2] < s.distancesToSecondCenters[site])
        {
            s.secondCenters[site] = site2;
            s.distancesToSecondCenters[site] = problem.distances[site][site2]; 
        }
    }
}