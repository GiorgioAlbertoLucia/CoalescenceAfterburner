#ifndef COALESCENCEENGINEHE3_H
#define COALESCENCEENGINEHE3_H

#include "CoalescenceEngine.h"
#include "WignerDensity.h"

class CoalescenceEngineHe3 : public CoalescenceEngine {
public:
    CoalescenceEngineHe3(Config& cfg) : CoalescenceEngine(cfg) { 
        fMass = fMassHe3;
        loadWignerDensity(cfg, fRadiusParameterHe3);
    }
    ~CoalescenceEngineHe3() {
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
    constexpr static int fA = 3; // Number of nucleons in He3
    constexpr static int fNRelativeCoords = fA - 1; // Number of relative Jacobi coordinates
    constexpr static int fNJacobiCoords = fA; // Total number of Jacobi coordinates (including CoM)
    constexpr static double fSA = 1/12.; // Spin-isospin factor for He3

    constexpr static double fMassHe3 = 2.80839160743; // [GeV/c^2] He3 mass
    constexpr static double fRadiusHe3 = 1.96; // [fm] He3 radius
    constexpr static double fRadiusParameterHe3 = fRadiusHe3 * std::sqrt(2.); // [fm] He3 radius parameter for Gaussian wavefunction

    WignerType fWignerType = WignerType::kSinglePair;
    WignerDensity* fWigner;
};

#endif // COALESCENCEENGINEHE3_H