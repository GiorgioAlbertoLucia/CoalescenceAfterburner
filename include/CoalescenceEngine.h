#ifndef COALESCENCEENGINE_H
#define COALESCENCEENGINE_H

#include "TDirectory.h"
#include "TH1D.h"

#include "BookKeeping.h"
#include "Config.h"
#include "Event.h"
#include "JacobiTransform.h"
#include "MomentumSampler.h"
#include "Particle.h"
#include "PositionSampler.h"

class CoalescenceEngine {
public:
    
    CoalescenceEngine(Config& config);
    virtual ~CoalescenceEngine() = default;
    void run(BookKeeping& bookKeeping);
    virtual void workerRun(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                           std::pair<float, float>& threadResult, std::vector<double>& threadYields, 
                           BookKeeping& bookKeeping) const = 0;
    void generateEvent(Event& event, PositionSampler& positionSampler, MomentumSampler& momentumSampler, TRandom3& random) const;
    virtual double runEvent(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform) const = 0;
    bool checkNucleusWithinRapidity(const std::vector<Particle>& nucleons, double yMin = -0.5, double yMax = 0.5) const;

protected:
    constexpr static double fMassNucleon = 0.9382720813; // [GeV/c^2] nucleon mass
    constexpr static double fYMax = 1.5; // half-width of the flat rapidity window for momentum sampling
    double fMass = 0.0; // mass of the nucleus (to be set in derived classes)

    Config          fConfig;
    double          fAverageNucleons = 0.0; // average number of nucleons per event
    TH1D*           fHPt{nullptr}; // input pT histogram (cloned from config file)
    BookKeeping     fBookKeeping;

};

#endif // COALESCENCEENGINE_H