#ifndef _ALGORITHM_H_
#define _ALGORITHM_H_ 

#include "structures.hpp"
#include <random>
#include <chrono>

class Algorithm
{
public:
    Algorithm(const Problem& problem, const std::mt19937& generator, int timeMax = 60 * 1000, bool verbose = false, double kstepCoef = 0.5, double kmaxCoef = 0.1, double ls3Prob = 0.52, double ls2Prob = 0.24)
        :problem(problem), generator(generator)
        , neighborhoodRuns({0,0,0}), neighborhoodImprovements({0,0,0})
        , timeMax(timeMax), verbose(verbose)
        , kstepCoef(kstepCoef), kmaxCoef(kmaxCoef), ls3Prob(ls3Prob), ls2Prob(ls2Prob)
    {}

    Result execute();

private:
    int getRunningTime() const;
    bool createGreedySolution(Solution& s);
    void localSearchSwap(Solution& s, bool& hasImproved);
    void localSearch2Out1In(Solution& s, bool& hasImproved);
    void localSearch1Out2In(Solution& s, bool& hasImproved);
    void localSearchSwapWithIndexRemembering(Solution& s, bool& hasImproved);
    void vnd(Solution& s, bool shuffle = false, bool doLS3 = true);
    void shake(Solution& s, int k);
    void shakeSimple(Solution& s, int k);
    void shake2(Solution& s, int k);
    void addAditionalSites(Solution& s, double additionalSiteProb);
    Result multistartVND();
    Result multistartVND2();
    Result VNS();

    DistanceAndCriticalCount getSwapEvaluation(int out, int in, const Solution& s);
    DistanceAndCriticalCount get2Out1InEvaluation(int out1, int out2, int in, const Solution& s);
    DistanceAndCriticalCount get1Out2InEvaluation(int out, int in1, int in2, const Solution& s);
    void performSwap(int outIdx, int inIdx, Solution& s);
    void perform2Out1In(int out1Idx, int out2Idx, int inIdx, Solution& s);
    void perform1Out2In(int outIdx, int in1Idx, int in2Idx, Solution& s);
    void calculateSecondCenter(int site, Solution& s);

    const Problem& problem;
    std::mt19937 generator;
    std::vector<int> neighborhoodRuns;
    std::vector<int> neighborhoodImprovements;
    std::chrono::_V2::steady_clock::time_point startTime;
    int timeMax;
    bool verbose;
    double kstepCoef;
    double kmaxCoef;
    double ls3Prob;
    double ls2Prob;
};

#endif