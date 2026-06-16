#ifndef PARTICLE_H
#define PARTICLE_H

#include "TLorentzVector.h"
#include "TVector3.h"

// PDG codes for the species we care about
namespace PDG {
    constexpr int kProton   = 2212;
    constexpr int kHe3      = 1000020030;
    constexpr int kLi4      = 1000030040;
}

struct Particle {

    int    pdg;         // PDG code
    TLorentzVector mom; // 4-momentum (px, py, pz, E) in GeV/c and GeV
    TVector3       pos; // emission position (x, y, z) in fm

    // --- convenience accessors ---

    // transverse momentum
    double pT() const { return mom.Pt(); }

    // rapidity
    double y() const { return mom.Rapidity(); }

    // transverse mass: mt = sqrt(E^2 - pz^2)
    // equivalently sqrt(m^2 + pT^2)
    double mT() const { return mom.Mt(); }

    // rest mass
    double mass() const { return mom.M(); }

    // --- constructors ---

    Particle() : pdg(0), mom(), pos() {}

    Particle(int pdg_, double px, double py, double pz, double E)
        : pdg(pdg_), mom(px, py, pz, E), pos(0., 0., 0.) {}
};

#endif // PARTICLE_H