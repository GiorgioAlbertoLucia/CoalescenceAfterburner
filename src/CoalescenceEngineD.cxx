#include "CoalescenceEngineD.h"
#include "CMFrameBooster.h"

void CoalescenceEngineD::workerRun(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
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

void CoalescenceEngineD::workerRunSource(long long evStart, long long evEnd, int threadId, TH1D* hPtNucleonClone,
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

double CoalescenceEngineD::runEvent(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, const WignerDensity* wigner) const {
    
    float result = 0.f; 
    if (event.nProtons() < 1 || event.nNeutrons() < 1) return result;
 
    std::vector<Particle> protons = event.protons;
    std::vector<Particle> neutrons = event.neutrons;
    bookKeeping.fHMultiplicityNucleons->Fill(protons.size());

    for (std::size_t iP1 = 0; iP1 < protons.size(); ++iP1) {

        if (true) {
            bookKeeping.fHPtNucleon->Fill(protons[iP1].mom.Pt());
            bookKeeping.fHYNucleon->Fill(protons[iP1].mom.Rapidity());
            bookKeeping.fHPhiNucleon->Fill(protons[iP1].mom.Phi());
            if (protons[iP1].mom.Rapidity() > -0.5 && protons[iP1].mom.Rapidity() < 0.5) {
                bookKeeping.fHPtProtonOneRapidityUnit->Fill(protons[iP1].mom.Pt());
            }
        }
        for (std::size_t iN1 = 0; iN1 < neutrons.size(); ++iN1) {
                    
            const Particle& p1 = protons[iP1];
            const Particle& n1 = neutrons[iN1];
            std::vector<Particle> nucleons = {p1, n1};

            if (!checkNucleusWithinRapidity(nucleons)) {
                continue; // skip if nucleus rapidity is outside [-0.5, 0.5]
            }

            std::vector<Particle> boostedNucleons = CMFrameBooster::boostParticles(nucleons);
            std::vector<Particle> jacobiNucleons = jacobiTransform.transform(boostedNucleons);
            
            // Evaluate A-body Wigner density
            double D = 1.;
            for (int j = 1; j < fA; ++j) {
                D *= fWigner->evaluate(jacobiNucleons[j].mom.Vect(), jacobiNucleons[j].pos);
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
 
    return result;
}

void CoalescenceEngineD::runEventAfterburner(Event& event, BookKeeping& bookKeeping, JacobiTransform& jacobiTransform, const WignerDensity* wigner, TRandom3& random) const {
    
    if (event.nProtons() < 1 || event.nNeutrons() < 1) return;
 
    std::vector<Particle> protons = event.protons;
    event.protonUsed.resize(protons.size(), false);
    std::vector<Particle> neutrons = event.neutrons;
    event.neutronUsed.resize(neutrons.size(), false);

    bookKeeping.fHMultiplicityNucleons->Fill(protons.size());
    
    for (std::size_t iP1 = 0; iP1 < protons.size(); ++iP1) {    
        
        if (true) {
            bookKeeping.fHPtNucleon->Fill(protons[iP1].mom.Pt());
            bookKeeping.fHYNucleon->Fill(protons[iP1].mom.Rapidity());
            bookKeeping.fHPhiNucleon->Fill(protons[iP1].mom.Phi());
            bookKeeping.fHPositionNucleons->Fill(protons[iP1].pos.Mag());
            if (protons[iP1].mom.Rapidity() > -0.5 && protons[iP1].mom.Rapidity() < 0.5) {
                bookKeeping.fHPtProtonOneRapidityUnit->Fill(protons[iP1].mom.Pt());
            }
        }
        if (event.protonUsed[iP1]) continue;

        for (std::size_t iN1 = 0; iN1 < neutrons.size(); ++iN1) {
            if (event.neutronUsed[iN1]) continue;
                    
            const Particle& p1 = protons[iP1];
            const Particle& n1 = neutrons[iN1];
            std::vector<Particle> nucleons = {p1, n1};


            
            std::vector<Particle> boostedNucleons = CMFrameBooster::boostParticles(nucleons);
            std::vector<Particle> jacobiNucleons = jacobiTransform.transform(boostedNucleons);
            Particle relativeParticle;
            relativeParticle.mom = boostedNucleons[0].mom - boostedNucleons[1].mom;
            relativeParticle.pos = 0.5 * (boostedNucleons[0].pos - boostedNucleons[1].pos);
            
            // Evaluate A-body Wigner density
            /// double D = 1.;
            /// double Dmax = 1.;
            /// for (int j = 1; j < fA; ++j) {
                ///     D *= wigner->evaluate(jacobiNucleons[j].mom.Vect(), jacobiNucleons[j].pos);
                ///     Dmax *= wigner->maxValue();
                /// }
                double D = wigner->evaluate(relativeParticle.mom.Vect(), relativeParticle.pos);
                double Dmax = wigner->maxValue();
                const double randomNumber = random.Rndm();
                if (randomNumber > D / Dmax) {
                    continue; // reject this combination based on Wigner density
                }
                
                // Mark nucleons as used
                event.protonUsed[iP1] = true;
                event.neutronUsed[iN1] = true;
                // Create the nucleus particle and add it to the event
                TLorentzVector pTot;
                for (const auto& p : nucleons)
                pTot += p.mom;
                Particle nucleus;
                nucleus.mom = pTot;
                nucleus.pos = 0.5 * (nucleons[0].pos + nucleons[1].pos);  // center of mass
                
                bookKeeping.fHPtNucleusBeforeYcut->Fill(pTot.Pt());

                if (!checkNucleusWithinRapidity(nucleons)) {
                    continue; // skip if nucleus rapidity is outside [-0.5, 0.5]
                }
                event.nuclei.push_back(nucleus);
                
                bookKeeping.fHPtNucleus->Fill(pTot.Pt());
                bookKeeping.fHYNucleus->Fill(pTot.Rapidity());
                bookKeeping.fHPhiNucleus->Fill(pTot.Phi());
                bookKeeping.fHPositionNucleus->Fill(nucleus.pos.Mag());
            }
    }
}
