#ifndef _ALGORITHM_H_
#define _ALGORITHM_H_ 

#include "structures.hpp"
#include <random>
#include <chrono>

class Algorithm
{
public:
    Algorithm(const Problem& problem, const std::mt19937& generator, int timeMax = 60 * 1000, bool verbose = false)
        :problem(problem), generator(generator)
        , neighborhoodRuns({0,0,0}), neighborhoodImprovements({0,0,0})
        , timeMax(timeMax), verbose(verbose)
    {}

    Result execute();

private:
    int getRunningTime() const;
    bool createGreedySolution(Solution& s);
    void localSearchSwap(Solution& s, bool& hasImproved);
    void localSearch2Out1In(Solution& s, bool& hasImproved);
    void localSearch1Out2In(Solution& s, bool& hasImproved);
    void localSearchSwapWithIndexRemembering(Solution& s, bool& hasImproved);
    void Algorithm::vnd(Solution& s, bool shuffle = false, bool doLS3 = true);
    void shake(Solution& s, int k);
    void shakeSimple(Solution& s, int k);
    void addAditionalSites(Solution& s, double additionalSiteProb);
    Result multistartVND();
    Result VNS();

    const Problem& problem;
    std::mt19937 generator;
    std::vector<int> neighborhoodRuns;
    std::vector<int> neighborhoodImprovements;
    std::chrono::_V2::steady_clock::time_point startTime;
    int timeMax;
    bool verbose;
};

#endif