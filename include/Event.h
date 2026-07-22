#ifndef EVENTA4_H
#define EVENTA4_H

#include "Particle.h"

#include <fstream>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// Event
//
// Plain container for one collision event.
// Holds only the species relevant for any A coalescence: protons and neutrons.
// ─────────────────────────────────────────────────────────────────────────────

struct Event {
    int                  id;       // event number as read from file
    std::vector<Particle> protons; // PDG 2212
    std::vector<Particle> neutrons;// PDG 2112

    void clear() {
        id = 0;
        protons.clear();
        neutrons.clear();
    }

    int nProtons() const { return static_cast<int>(protons.size()); }
    int nNeutrons() const { return static_cast<int>(neutrons.size()); }
};

#endif // EVENTA4_H