#ifndef EVENTREADER_H
#define EVENTREADER_H

#include "Particle.h"

#include <fstream>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// Event
//
// Plain container for one collision event.
// Holds only the species relevant for Li4 coalescence: He3 and protons.
// All other particles in the FIST output are skipped during parsing.
// ─────────────────────────────────────────────────────────────────────────────

struct Event {
    int                  id;       // event number as read from file
    std::vector<Particle> protons; // PDG 2212
    std::vector<Particle> he3s;    // PDG 1000020030

    void clear() {
        id = 0;
        protons.clear();
        he3s.clear();
    }

    int nProtons() const { return static_cast<int>(protons.size()); }
    int nHe3()     const { return static_cast<int>(he3s.size());    }
};


// ─────────────────────────────────────────────────────────────────────────────
// EventReader
//
// Reads the FIST-sampler ASCII output format event by event.
//
// Expected file structure:
//
//   Event 1
//           pdgid    p0[GeV/c2]    px[GeV/c]    py[GeV/c]    pz[GeV/c]
//   2212    0.938    0.312    -0.124    0.087
//   1000020030    2.809    -0.201    0.443    -0.312
//   ...
//   Event 2
//   ...
//
// Notes:
//   - The header line containing column labels is skipped automatically.
//   - Particles with PDG codes other than proton or He3 are discarded.
//   - Antiparticles (PDG -2212, -1000020030) are also discarded by default
//     but can be included via setReadAntiparticles(true).
//   - p0 (energy) is read but recomputed from (px,py,pz,m) for consistency.
//
// Usage:
//   EventReader reader("fist_output.dat");
//   Event ev;
//   while (reader.nextEvent(ev)) {
//       // process ev.protons and ev.he3s
//   }
// ─────────────────────────────────────────────────────────────────────────────

class EventReader {
public:

    explicit EventReader(const std::string& filename);
    ~EventReader();

    // Reads the next event into ev. Returns true if successful, false on EOF.
    bool nextEvent(Event& ev);

    // Rewind to the beginning of the file
    void rewind();

    // Total events read so far
    int eventsRead() const { return fEventsRead; }

    // Include antiparticles (PDG -2212, -1000020030) in the output.
    // When true, antiprotons go into ev.protons and anti-He3 into ev.he3s.
    // Default: false.
    void setReadAntiparticles(bool flag) { fReadAnti = flag; }

    bool isOpen() const { return fFile.is_open(); }

private:

    // Parse a particle line and append to ev if it is a relevant species.
    // Returns false if the line could not be parsed.
    bool parseLine(const std::string& line, Event& ev) const;

    // Check whether a PDG code (absolute value) is one we care about
    static bool isRelevant(int absPdg) {
        return absPdg == PDG::kProton || absPdg == PDG::kHe3;
    }

    std::string   fFilename;
    std::ifstream fFile;
    int           fEventsRead;
    bool          fReadAnti;

    // Buffer for the line already consumed to detect the next "Event N" header
    std::string   fPendingLine;
    bool          fHasPending;
};

#endif // EVENTREADER_H