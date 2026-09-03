#ifndef EVENTA4_H
#define EVENTA4_H

#include "Particle.h"

#include <fstream>
#include <string>
#include <vector>

/**
 * Plain container for one collision event.
 * Holds only the species relevant for any A coalescence: protons and neutrons.
 */
struct Event {
    int                  id;       // event number as read from file
    std::vector<Particle> protons; // PDG 2212
    std::vector<Particle> neutrons;// PDG 2112
    std::vector<Particle> nuclei;  // nucleus produced with the coalescence, for afterburner output
    std::vector<bool> protonUsed; // flags to indicate if a proton has been used in coalescence
    std::vector<bool> neutronUsed;// flags to indicate if a neutron has been used in coalescence

    void clear() {
        id = 0;
        protons.clear();
        neutrons.clear();
        nuclei.clear();
    }

    int nProtons() const { return static_cast<int>(protons.size()); }
    int nNeutrons() const { return static_cast<int>(neutrons.size()); }
};

#endif // EVENTA4_H