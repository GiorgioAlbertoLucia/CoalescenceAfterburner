#include "CoalescenceEngineHe3.h"
#include "CMFrameBooster.h"

void CoalescenceEngineHe3::workerRun(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                                     std::pair<float, float>& threadResult, std::vector<double>& threadYields,
                                     BookKeeping& bookKeeping) const {

    bookKeeping.init(threadId);
    PositionSampler positionSampler(fConfig.sourceRadius, fConfig.randomSeed + 1 * threadId);
    MomentumSampler momentumSampler(fMassNucleon, hPtNucleonClone, fYMax, fConfig.randomSeed + 2 * threadId);
    TRandom3 random(fConfig.randomSeed + 3 * threadId);
    JacobiTransform jacobiTransform(fA);
    float yield = 0.f;

    const long long nLocal   = evEnd - evStart;
    const long long printEvery = std::max(1LL, nLocal / 5);
 
    for (long long iEv = evStart; iEv < evEnd; ++iEv) {
        if ((iEv - evStart) % printEvery == 0)
            std::cout << " [Thread " << threadId << "] Event "
                      << (iEv - evStart) << " / " << nLocal
                      << "  (" << (5 * (iEv - evStart) / nLocal) << "%)\n";

        Event event;
        generateEvent(event, positionSampler, momentumSampler, random);

        const float weight = runEvent(event, bookKeeping, jacobiTransform);
        yield += weight;
        threadYields.push_back(weight);
    }

    threadResult = std::make_pair(yield, static_cast<float>(nLocal));
    std::cout << " [Thread " << threadId << "] Finished " << nLocal
              << " events, yield = " << yield << "\n";
}

double CoalescenceEngineHe3::runEvent(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform) const {
    
    float result = 0.f; 
    if (event.nProtons() < 2 || event.nNeutrons() < 1) return result;
 
    std::vector<Particle> protons = event.protons;
    std::vector<Particle> neutrons = event.neutrons;
 
    // Loop over all 3 nucleon combinations (2p + 1n) and apply coalescence.
    // All particles can coalesce any time (no "used" tracking), so we just sum the weights.
    // Always use different protons (no self-pairs)
    for (std::size_t iP1 = 0; iP1 < protons.size(); ++iP1) {
        for (std::size_t iP2 = iP1 + 1; iP2 < protons.size(); ++iP2) {
            for (std::size_t iN1 = 0; iN1 < neutrons.size(); ++iN1) {
                    
                const Particle& p1 = protons[iP1];
                const Particle& p2 = protons[iP2];
                const Particle& n1 = neutrons[iN1];
                std::vector<Particle> nucleons = {p1, p2, n1};

                if (true) {
                    for (const auto& p : nucleons) {
                        bookKeeping.fHPtNucleon->Fill(p.mom.Pt());
                        bookKeeping.fHYNucleon->Fill(p.mom.Rapidity());
                        bookKeeping.fHPhiNucleon->Fill(p.mom.Phi());
                    }
                }


                if (!checkNucleusWithinRapidity(nucleons)) {
                    continue; // skip if nucleus rapidity is outside [-0.5, 0.5]
                }

                std::vector<Particle> boostedNucleons = CMFrameBooster::boostParticles(nucleons);
                std::vector<Particle> jacobiNucleons = jacobiTransform.transform(boostedNucleons);
                
                // Evaluate A-body Wigner density
                double D = 1.;
                for (int j = 1; j < fA; ++j) {
                    D *= fWignerSinglePair->evaluate(jacobiNucleons[j].mom.Vect(), jacobiNucleons[j].pos);
                }
            
                result += static_cast<float>(fSA * D);
                if (true) {
                    TLorentzVector pTot;
                    for (const auto& p : nucleons)
                        pTot += p.mom;
                    bookKeeping.fHPtNucleus->Fill(pTot.Pt(), result);
                    bookKeeping.fHYNucleus->Fill(pTot.Rapidity(), result);
                    bookKeeping.fHPhiNucleus->Fill(pTot.Phi(), result);
                }
            }
        }
    }
 
    return result;
}
