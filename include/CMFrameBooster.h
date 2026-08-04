#ifndef CMFRAMEBOOSTER_H
#define CMFRAMEBOOSTER_H

#include "Particle.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <vector>
#include <stdexcept>

/**
 * @brief Boosts a set of A particles from the lab frame into their centre-of-mass
 * (CM) frame, returning the boosted 3-momenta and 3-positions separately.
 *
 * Momentum boost:
 *   The CM frame beta vector is derived from the total 4-momentum:
 *     beta = P_tot.BoostVector()
 *   Each particle 4-momentum is then boosted by -beta.
 *   The output 3-momenta are the spatial components of the boosted 4-momenta.
 *
 * Position boost:
 *   Positions are treated as 4-vectors with t=0 (equal-time approximation,
 *   consistent with the paper Sec. 2 and CoalescenceEngine::boostToPRF).
 *   Each position 4-vector (x, y, z, 0) is boosted by -beta.
 *   The output 3-positions are the spatial components of the boosted vectors.
 *
 * Usage:
 *   auto [momCM, posCM] = CMFrameBooster::boost(particles);
 *   // then pass momCM and posCM to JacobiTransform::relative()
 */
class CMFrameBooster {
public:

    static std::vector<Particle> boostParticles(const std::vector<Particle>& particles) {
        if (particles.empty())
            throw std::invalid_argument("CMFrameBooster::boost: empty particle list");

        std::vector<Particle> boostedParticles;

        TLorentzVector pTot;
        for (const auto& p : particles) {
            pTot += p.mom;
        }
        const TVector3 beta = pTot.BoostVector();

        for (const auto& p : particles) {

            TLorentzVector boostedMom = p.mom;
            boostedMom.Boost(-beta);

            TLorentzVector boostedPos(p.pos[0], p.pos[1], p.pos[2], 0.);
            boostedPos.Boost(-beta);

            Particle boostedParticle;
            boostedParticle.pdg = p.pdg;
            boostedParticle.mom = boostedMom;
            boostedParticle.pos = TVector3(boostedPos.X(), boostedPos.Y(), boostedPos.Z());
            boostedParticles.push_back(boostedParticle);
            
            //boostedParticles.emplace_back(p.pdg, boostedMom, TVector3(boostedPos.X(), boostedPos.Y(), boostedPos.Z()));
        }

        return boostedParticles;
    }

};

#endif // CMFRAMEBOOSTER_H