#ifndef COALESCENCEENGINE_D_H
#define COALESCENCEENGINE_D_H

#include "CoalescenceEngine.h"
#include "WignerDensity.h"

class CoalescenceEngineD : public CoalescenceEngine {
public:
    CoalescenceEngineD(Config& cfg) : CoalescenceEngine(cfg) { 
        fMass = fMassD;
        loadWignerDensity(cfg, fRadiusParameterD);
    }
    ~CoalescenceEngineD() {
        if (fWigner) delete fWigner;
    }
    void workerRun(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                   std::pair<float, float>& threadResult, std::vector<double>& threadYields, 
                   const WignerDensity* wigner, BookKeeping& bookKeeping) const override;
    void workerRunSource(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                        const WignerDensity* wigner, BookKeeping& bookKeeping) const override;
    double runEvent(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, const WignerDensity* wigner) const override;
    void runEventAfterburner(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, 
                             const WignerDensity* wigner, TRandom3& random) const override;

private:
    constexpr static int fA = 2; // Number of nucleons in D
    constexpr static int fNRelativeCoords = fA - 1; // Number of relative Jacobi coordinates
    constexpr static int fNJacobiCoords = fA; // Total number of Jacobi coordinates (including CoM)
    constexpr static double fSA = 3/8.; // Spin-isospin factor for D

    constexpr static double fMassD = 1.87561294257; // [GeV/c^2] D mass
    constexpr static double fRadiusD = 2.13; // [fm] D radius
    constexpr static double fRadiusParameterD = fRadiusD * std::sqrt(8./3.); // [fm] D radius parameter for Gaussian wavefunction

    WignerType fWignerType = WignerType::kSinglePair;
    WignerDensity* fWigner;
};

#endif // COALESCENCEENGINE_D_H