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

    int    pdg;         
    TLorentzVector mom; // 4-momentum (px, py, pz, E) in GeV/c and GeV
    TVector3       pos; // emission position (x, y, z) in fm

    double pT() const { return mom.Pt(); }
    double y() const { return mom.Rapidity(); }
    double mT() const { return mom.Mt(); }
    double mass() const { return mom.M(); }

    Particle() : pdg(0), mom(), pos() {}
    Particle(int pdg_, double px, double py, double pz, double E)
        : pdg(pdg_), mom(px, py, pz, E), pos(0., 0., 0.) {}
    Particle(int pdg_, const TLorentzVector& mom_, const TVector3& pos_)
        : pdg(pdg_), mom(mom_), pos(pos_) {}
};

#endif // PARTICLE_H