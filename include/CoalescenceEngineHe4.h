#ifndef COALESCENCEENGINEHE4_H
#define COALESCENCEENGINEHE4_H

#include "CoalescenceEngine.h"
#include "WignerDensity.h"

class CoalescenceEngineHe4 : public CoalescenceEngine {
public:
    CoalescenceEngineHe4(Config& cfg) : CoalescenceEngine(cfg) { 
        fMass = fMassHe4;
        loadWignerDensity(cfg, fRadiusParameterHe4);
    }
    ~CoalescenceEngineHe4() {
        if (fWigner) delete fWigner;
    }
    void workerRun(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                   std::pair<float, float>& threadResult, std::vector<double>& threadYields, 
                   const WignerDensity* wigner, BookKeeping& bookKeeping) const override;
    void workerRunSource(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                        const WignerDensity* wigner, BookKeeping& bookKeeping) const override { throw std::runtime_error("workerRunSource is not implemented for CoalescenceEngineHe4"); };
    double runEvent(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, const WignerDensity* wigner) const override;
    void runEventAfterburner(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, const WignerDensity* wigner, TRandom3& random) const override {
        throw std::runtime_error("runEventAfterburner is not implemented for CoalescenceEngineHe4");
    }

private:
    constexpr static int fA = 4; // Number of nucleons in He4
    constexpr static int fNRelativeCoords = fA - 1; // Number of relative Jacobi coordinates
    constexpr static int fNJacobiCoords = fA; // Total number of Jacobi coordinates (including CoM)
    constexpr static double fSA = 1/96.; // Spin-isospin factor for He4

    constexpr static double fMassHe4 = 3.7273794066; // [GeV/c^2] He4 mass
    constexpr static double fRadiusHe4 = 1.67; // [fm] He4 radius
    constexpr static double fRadiusParameterHe4 = fRadiusHe4 * std::sqrt(2. * 2./3. * 4./3.); // [fm] He4 radius parameter for Gaussian wavefunction
    
    WignerType fWignerType = WignerType::kSinglePair;
    WignerDensity* fWigner;
};

#endif // COALESCENCEENGINEHE4_H