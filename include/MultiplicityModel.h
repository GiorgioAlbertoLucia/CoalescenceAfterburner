#ifndef MULTIPLICITYMODEL_H
#define MULTIPLICITYMODEL_H

#include "TRandom3.h"

#include <cmath>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// MultiplicityModel
//
// Determines the number of protons and He3s to generate per event,
// following the ToMCCA philosophy: the particle composition is driven by
// the ratio of the integrated yields of each species.
//
// The number of particles per event is drawn from one of three modes:
//
//   Fixed      : exactly <N> particles every event (useful for tests)
//   Poissonian : N ~ Poisson(<N>)  (simple, fast)
//   External   : caller provides N directly (e.g. from a multiplicity file)
//
// The mean number of protons and He3s per event is provided directly as
// <Np> and <NHe3>, computed by the caller from the integral of the input
// pT spectra scaled by the rapidity window and the number of events.
//
// In the ToMCCA spirit, p and He3 multiplicities are treated independently
// (uncorrelated emission), which is appropriate for a two-species system
// where one species is already a composite nucleus.
// ─────────────────────────────────────────────────────────────────────────────

enum class MultiplicityMode {
    Fixed,
    Poissonian
};

class MultiplicityModel {
public:

    // meanNp   : mean number of protons per event
    // meanNHe3 : mean number of He3 per event
    // mode     : Fixed or Poissonian
    // seed     : TRandom3 seed
    MultiplicityModel(double           meanNp,
                      double           meanNHe3,
                      MultiplicityMode mode = MultiplicityMode::Poissonian,
                      unsigned int     seed = 43)
        : fMeanNp(meanNp)
        , fMeanNHe3(meanNHe3)
        , fMode(mode)
        , fRng(seed)
    {
        if (meanNp < 0.)
            throw std::invalid_argument(
                "MultiplicityModel: meanNp must be non-negative");
        if (meanNHe3 < 0.)
            throw std::invalid_argument(
                "MultiplicityModel: meanNHe3 must be non-negative");
    }

    // Draw number of protons for one event
    int drawNProtons() const  { return draw(fMeanNp);   }

    // Draw number of He3s for one event
    int drawNHe3() const      { return draw(fMeanNHe3); }

    void setSeed(unsigned int seed) { fRng.SetSeed(seed); }

    double meanNProtons() const { return fMeanNp;   }
    double meanNHe3()     const { return fMeanNHe3; }

private:
    int draw(double mean) const {
        if (mean <= 0.) return 0;
        switch (fMode) {
            case MultiplicityMode::Fixed:
                return static_cast<int>(std::round(mean));
            case MultiplicityMode::Poissonian:
                return static_cast<int>(fRng.Poisson(mean));
        }
        return 0;
    }

    double           fMeanNp;
    double           fMeanNHe3;
    MultiplicityMode fMode;
    mutable TRandom3 fRng;
};

#endif // MULTIPLICITYMODEL_H
