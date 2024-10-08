#ifndef _ALGORITHM_H_
#define _ALGORITHM_H_ 

#include "structures.hpp"
#include <random>
#include <chrono>

//     kstep kmax ls3prob ls2prob
// GVNS: 0.5 0.1 0.52 0.24 
// GVNS-LS3 : 0.3 0.1 0.25 0.40
// 
// GVNS: Algorithm(const Problem& problem, const std::mt19937& generator, int timeMax = 60 * 1000, bool verbose = false, bool doLS3 = false, double kstepCoef = 0.5, double kmaxCoef = 0.1, double ls3Prob = 0.52, double ls2Prob = 0.24)
// GVNS-LS3: Algorithm(const Problem& problem, const std::mt19937& generator, int timeMax = 60 * 1000, bool verbose = false, bool doLS3 = true, double kstepCoef = 0.3, double kmaxCoef = 0.1, double ls3Prob = 0.25, double ls2Prob = 0.40)

class Algorithm
{
public:  
    Algorithm(const Problem& problem, const std::mt19937& generator, int timeMax = 60 * 1000, bool verbose = false, bool doLS3 = false, double kstepCoef = 0.5, double kmaxCoef = 0.1, double ls3Prob = 0.33/*0.52*/, double ls2Prob = 0.27/*0.24*/)
        :problem(problem), generator(generator)
        , neighborhoodRuns({0,0,0}), neighborhoodImprovements({0,0,0})
        , timeMax(timeMax), verbose(verbose), doLS3(doLS3)
        , kstepCoef(kstepCoef), kmaxCoef(kmaxCoef), ls3Prob(ls3Prob), ls2Prob(ls2Prob)
    {}

    Result execute();

private:
    int getRunningTime() const;
    bool createGreedySolution(Solution& s);
    void localSearchSwap(Solution& s, bool& hasImproved);
    void localSearch2Out1In(Solution& s, bool& hasImproved);
    void localSearch1Out2In(Solution& s, bool& hasImproved);
    void vnd(Solution& s, bool shuffle = false, bool doLS3 = true);
    void shake(Solution& s, int k);
    DistanceAndCriticalCount getSwapEvaluation(int out, int in, const Solution& s);
    DistanceAndCriticalCount get2Out1InEvaluation(int out1, int out2, int in, const Solution& s);
    DistanceAndCriticalCount get1Out2InEvaluation(int out, int in1, int in2, const Solution& s);
    void performSwap(int outIdx, int inIdx, Solution& s);
    void perform2Out1In(int out1Idx, int out2Idx, int inIdx, Solution& s);
    void perform1Out2In(int outIdx, int in1Idx, int in2Idx, Solution& s);
    void calculateSecondCenter(int site, Solution& s);
    void shakeSwap(Solution& s);
    void shake2Out1In(Solution& s);
    void shake1Out2In(Solution& s);
    Result multistartVND();
    Result VNS();

    const Problem& problem;
    std::mt19937 generator;
    std::vector<int> neighborhoodRuns;
    std::vector<int> neighborhoodImprovements;
    std::chrono::_V2::steady_clock::time_point startTime;
    int timeMax;
    bool verbose;
    bool doLS3;
    double kstepCoef;
    double kmaxCoef;
    double ls3Prob;
    double ls2Prob;
};

#endif