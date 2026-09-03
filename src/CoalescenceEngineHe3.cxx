#include "CoalescenceEngineHe3.h"
#include "CMFrameBooster.h"

void CoalescenceEngineHe3::workerRun(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                                     std::pair<float, float>& threadResult, std::vector<double>& threadYields,
                                     const WignerDensity* wigner, BookKeeping& bookKeeping) const {

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

        const float weight = runEvent(event, bookKeeping, jacobiTransform, wigner);
        yield += weight;
        threadYields.push_back(weight);
    }

    threadResult = std::make_pair(yield, static_cast<float>(nLocal));
    std::cout << " [Thread " << threadId << "] Finished " << nLocal
              << " events, yield = " << yield << "\n";
}

void CoalescenceEngineHe3::workerRunSource(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
                                     const WignerDensity* wigner, BookKeeping& bookKeeping) const {

    bookKeeping.init(threadId);
    PositionSampler positionSampler(fConfig.sourceRadius, fConfig.randomSeed + 1 * threadId);
    MomentumSampler momentumSampler(fMassNucleon, hPtNucleonClone, fYMax, fConfig.randomSeed + 2 * threadId);
    TRandom3 random(fConfig.randomSeed + 3 * threadId);
    JacobiTransform jacobiTransform(fA);

    const long long nLocal   = evEnd - evStart;
    const long long printEvery = std::max(1LL, nLocal / 5);
 
    for (long long iEv = evStart; iEv < evEnd; ++iEv) {
        if ((iEv - evStart) % printEvery == 0)
            std::cout << " [Thread " << threadId << "] Event "
                      << (iEv - evStart) << " / " << nLocal
                      << "  (" << (5 * (iEv - evStart) / nLocal) << "%)\n";

        Event event;
        generateEvent(event, positionSampler, momentumSampler, random);
        runEventAfterburner(event, bookKeeping, jacobiTransform, wigner, random);

        for (int iP = 0; iP < event.nProtons(); ++iP) {
            if (event.protonUsed[iP]) {
                continue;
            }
            for (const auto& nucleus : event.nuclei) {
                Particle protonCopy = event.protons[iP];
                Particle nucleusCopy = nucleus;
                std::vector<Particle> pair = {protonCopy, nucleusCopy};
                std::vector<Particle> boostedPair = CMFrameBooster::boostParticles(pair);

                const double r = (boostedPair[0].pos - boostedPair[1].pos).Mag();
                bookKeeping.fHRelativePosition->Fill(r);
            }
        }
    }

    std::cout << " [Thread " << threadId << "] Finished " << nLocal
              << " events\n";
}

double CoalescenceEngineHe3::runEvent(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, const WignerDensity* wigner) const {
    
    float result = 0.f; 
    if (event.nProtons() < 2 || event.nNeutrons() < 1) return result;
 
    std::vector<Particle> protons = event.protons;
    std::vector<Particle> neutrons = event.neutrons;

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
                if (fWignerType == WignerType::kSinglePair) {
                    for (int j = 1; j < fA; ++j) {
                        D *= wigner->evaluate(jacobiNucleons[j].mom.Vect(), jacobiNucleons[j].pos);
                    }
                } else if (fWignerType == WignerType::kFull) {
                    throw std::runtime_error("Full Wigner density is not implemented for CoalescenceEngineHe3");
                } else {
                    throw std::runtime_error("Unknown WignerType in CoalescenceEngineHe3::runEvent");
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

void CoalescenceEngineHe3::runEventAfterburner(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, const WignerDensity* wigner, TRandom3& random) const {
    
    if (event.nProtons() < 1 || event.nNeutrons() < 1) return;
 
    std::vector<Particle> protons = event.protons;
    event.protonUsed.resize(protons.size(), false);
    std::vector<Particle> neutrons = event.neutrons;
    event.neutronUsed.resize(neutrons.size(), false);

    for (std::size_t iP1 = 0; iP1 < protons.size(); ++iP1) {
        if (event.protonUsed[iP1]) continue;

        for (std::size_t iP2 = iP1 + 1; iP2 < protons.size(); ++iP2) {
            if (event.protonUsed[iP2]) continue;

            for (std::size_t iN1 = 0; iN1 < neutrons.size(); ++iN1) {
                if (event.neutronUsed[iN1]) continue;
                    
                const Particle& p1 = protons[iP1];
                const Particle& p2 = protons[iP2];
                const Particle& n1 = neutrons[iN1];
                std::vector<Particle> nucleons = {p1, p2, n1};

                if (true) {
                    for (const auto& p : nucleons) {
                        bookKeeping.fHPtNucleon->Fill(p.mom.Pt());
                        bookKeeping.fHYNucleon->Fill(p.mom.Rapidity());
                        bookKeeping.fHPhiNucleon->Fill(p.mom.Phi());
                        bookKeeping.fHPositionNucleons->Fill(p.pos.Mag());
                    }
                }

                if (!checkNucleusWithinRapidity(nucleons)) {
                    continue; // skip if nucleus rapidity is outside [-0.5, 0.5]
                }

                std::vector<Particle> boostedNucleons = CMFrameBooster::boostParticles(nucleons);
                std::vector<Particle> jacobiNucleons = jacobiTransform.transform(boostedNucleons);
                
                // Evaluate A-body Wigner density
                double D = 1.;
                double Dmax = 1.;
                if (fWignerType == WignerType::kSinglePair) {
                    for (int j = 1; j < fA; ++j) {
                        D *= wigner->evaluate(jacobiNucleons[j].mom.Vect(), jacobiNucleons[j].pos);
                        Dmax *= wigner->maxValue();
                    }
                } else if (fWignerType == WignerType::kFull) {
                    // adjust to support full Wigner density
                    // TODO: Implement full Wigner density evaluation
                    D = 1.;
                    Dmax = 1.;
                } else {
                    throw std::runtime_error("Unknown WignerType in CoalescenceEngineHe3::runEvent");
                }

                // Mark nucleons as used
                event.protonUsed[iP1] = true;
                event.protonUsed[iP2] = true;
                event.neutronUsed[iN1] = true;
                // Create the nucleus particle and add it to the event
                TLorentzVector pTot;
                for (const auto& p : nucleons)
                    pTot += p.mom;
                Particle nucleus;
                nucleus.mom = pTot;
                nucleus.pos = 0.5 * (nucleons[0].pos + nucleons[1].pos);  // center of mass
                event.nuclei.push_back(nucleus);
                
                bookKeeping.fHPtNucleus->Fill(pTot.Pt());
                bookKeeping.fHYNucleus->Fill(pTot.Rapidity());
                bookKeeping.fHPhiNucleus->Fill(pTot.Phi());
                bookKeeping.fHPositionNucleus->Fill(nucleus.pos.Mag());
            }
        }
    }
}

