#ifndef CMFRAMEBOOSTER_H
#define CMFRAMEBOOSTER_H

#include "Particle.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <vector>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// CMFrameBooster
//
// Boosts a set of A particles from the lab frame into their centre-of-mass
// (CM) frame, returning the boosted 3-momenta and 3-positions separately.
//
// Momentum boost:
//   The CM frame beta vector is derived from the total 4-momentum:
//     beta = P_tot.BoostVector()
//   Each particle 4-momentum is then boosted by -beta.
//   The output 3-momenta are the spatial components of the boosted 4-momenta.
//
// Position boost:
//   Positions are treated as 4-vectors with t=0 (equal-time approximation,
//   consistent with the paper Sec. 2 and CoalescenceEngine::boostToPRF).
//   Each position 4-vector (x, y, z, 0) is boosted by -beta.
//   The output 3-positions are the spatial components of the boosted vectors.
//
// Usage:
//   auto [momCM, posCM] = CMFrameBooster::boost(particles);
//   // then pass momCM and posCM to JacobiTransform::relative()
// ─────────────────────────────────────────────────────────────────────────────

class CMFrameBooster {
public:

    struct Result {
        std::vector<TVector3> momenta;   // 3-momenta in CM frame [GeV/c]
        std::vector<TVector3> positions; // 3-positions in CM frame [fm]
        TVector3              beta;      // boost vector used (lab -> CM)
    };

    // Boost all particles to the CM frame.
    // particles must have mom and pos filled.
    static Result boost(const std::vector<Particle>& particles) {
        if (particles.empty())
            throw std::invalid_argument("CMFrameBooster::boost: empty particle list");

        // ── Total 4-momentum ─────────────────────────────────────────────────
        TLorentzVector pTot;
        for (const auto& p : particles)
            pTot += p.mom;

        const TVector3 beta = pTot.BoostVector();

        // ── Boost each particle ───────────────────────────────────────────────
        Result result;
        result.beta = beta;
        result.momenta.reserve(particles.size());
        result.positions.reserve(particles.size());

        for (const auto& p : particles) {
            // Boost 4-momentum
            TLorentzVector mom = p.mom;
            mom.Boost(-beta);
            result.momenta.push_back(mom.Vect());

            // Boost position as a 4-vector with t=0 (equal-time approximation)
            TLorentzVector pos(p.pos, 0.);
            pos.Boost(-beta);
            result.positions.push_back(pos.Vect());
        }

        return result;
    }

    static std::vector<Particle> boostParticles(const std::vector<Particle>& particles) {
        if (particles.empty())
            throw std::invalid_argument("CMFrameBooster::boost: empty particle list");

        // ── Total 4-momentum ─────────────────────────────────────────────────
        TLorentzVector pTot;
        TVector3 posCM(0., 0., 0.);
        for (const auto& p : particles) {
            pTot += p.mom;
            posCM += p.pos;
        }

        const TVector3 beta = pTot.BoostVector();

        // ── Boost each particle ───────────────────────────────────────────────
        Result result;
        result.beta = beta;
        result.momenta.reserve(particles.size());
        result.positions.reserve(particles.size());

        for (const auto& p : particles) {
            // Boost 4-momentum
            TLorentzVector mom = p.mom;
            mom.Boost(-beta);
            result.momenta.push_back(mom.Vect());

            // Boost position as a 4-vector with t=0 (equal-time approximation)
            TLorentzVector pos(p.pos[0], p.pos[1], p.pos[2], 0.);
            pos.Boost(-beta);
            result.positions.push_back(pos.Vect());
        }

        // Build boosted particles
        std::vector<Particle> boostedParticles;
        boostedParticles.reserve(particles.size());
        for (std::size_t i = 0; i < particles.size(); ++i) {
            Particle p = particles[i];
            p.mom.SetVect(result.momenta[i]);
            const double mass = p.mom.M();
            p.mom.SetE(std::sqrt(result.momenta[i].Mag2() + mass * mass));
            p.pos = result.positions[i];
            boostedParticles.push_back(p);
        }
        return boostedParticles;
    }

    // Convenience overload: boost bare 3-vectors given a mass for each particle.
    // Useful when you have momenta and positions separately rather than Particles.
    static Result boost(const std::vector<TVector3>& momenta_lab,
                        const std::vector<TVector3>& positions_lab,
                        const std::vector<double>&   masses_GeV)
    {
        if (momenta_lab.size() != positions_lab.size())
            throw std::invalid_argument(
                "CMFrameBooster::boost: momenta and positions size mismatch");

        // Build Particle list then delegate
        std::vector<Particle> particles;
        particles.reserve(momenta_lab.size());
        for (std::size_t i = 0; i < momenta_lab.size(); ++i) {
            const double px = momenta_lab[i].X();
            const double py = momenta_lab[i].Y();
            const double pz = momenta_lab[i].Z();
            const double E  = std::sqrt(px*px + py*py + pz*pz
                                        + masses_GeV[i]*masses_GeV[i]);
            Particle p;
            p.mom.SetPxPyPzE(px, py, pz, E);
            p.pos = positions_lab[i];
            particles.push_back(p);
        }
        return boost(particles);
    }
};

#endif // CMFRAMEBOOSTER_H