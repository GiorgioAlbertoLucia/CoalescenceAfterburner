#include "EventReader.h"

#include "TLorentzVector.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <cstdlib> // std::abs for int

// PDG masses [GeV/c^2] used to recompute E from momentum
namespace {
    constexpr double kMassProton = 0.938272;
    constexpr double kMassHe3   = 2.808391;

    double massForPdg(int pdg) {
        const int absPdg = std::abs(pdg);
        if (absPdg == PDG::kProton) return kMassProton;
        if (absPdg == PDG::kHe3)   return kMassHe3;
        return 0.;
    }

    // Returns true if the string starts with "Event" (case-sensitive)
    bool isEventHeader(const std::string& line) {
        return line.rfind("Event", 0) == 0;
    }

    // Returns true if the line looks like the column-label header
    // (contains "pdgid" as a substring)
    bool isColumnHeader(const std::string& line) {
        return line.find("pdgid") != std::string::npos;
    }
}

// ─────────────────────────────────────────────────────────────────────────────

EventReader::EventReader(const std::string& filename)
    : fFilename(filename)
    , fEventsRead(0)
    , fReadAnti(true) // default to reading antiparticles, can be changed via setReadAntiparticles()
    , fHasPending(false)
{
    fFile.open(filename);
    if (!fFile.is_open())
        throw std::runtime_error("EventReader: cannot open file: " + filename);
}

EventReader::~EventReader() {
    if (fFile.is_open()) fFile.close();
}

void EventReader::rewind() {
    fFile.clear();
    fFile.seekg(0, std::ios::beg);
    fEventsRead  = 0;
    fHasPending  = false;
    fPendingLine = "";
}

// ─────────────────────────────────────────────────────────────────────────────

bool EventReader::nextEvent(Event& ev) {
    ev.clear();

    if (!fFile.is_open() || fFile.eof()) return false;

    std::string line;

    // ── 1. Find the "Event N" header ────────────────────────────────────────
    // On the first call fHasPending is false and we scan forward.
    // On subsequent calls the previous iteration already consumed one line
    // past the event boundary and stored it in fPendingLine.

    if (fHasPending) {
        line = fPendingLine;
        fHasPending = false;
    } else {
        // Scan until we find an event header or hit EOF
        while (std::getline(fFile, line)) {
            if (isEventHeader(line)) break;
        }
        if (fFile.eof() && !isEventHeader(line)) return false;
    }

    // Parse the event number from "Event N"
    {
        std::istringstream ss(line);
        std::string tag;
        ss >> tag >> ev.id;
    }

    // ── 2. Skip the column-label header line ────────────────────────────────
    if (!std::getline(fFile, line)) return false; // EOF right after header
    // Tolerate files where the column line is absent (e.g. blank line)
    if (!isColumnHeader(line) && !line.empty()) {
        // It is already a particle line or next event — treat as pending
        if (isEventHeader(line)) {
            fPendingLine = line;
            fHasPending  = true;
            ++fEventsRead;
            return true; // empty event
        }
        parseLine(line, ev); // it was a particle line
    }

    // ── 3. Read particle lines until next "Event" header or EOF ─────────────
    while (std::getline(fFile, line)) {
        if (line.empty()) continue;

        if (isEventHeader(line)) {
            // Store for next call and stop
            fPendingLine = line;
            fHasPending  = true;
            break;
        }

        if (isColumnHeader(line)) continue; // skip stray label lines

        parseLine(line, ev);
    }

    ++fEventsRead;
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────

bool EventReader::parseLine(const std::string& line, Event& ev) const {
    std::istringstream ss(line);

    int    pdg;
    double p0, px, py, pz;

    if (!(ss >> pdg >> p0 >> px >> py >> pz)) return false;

    const int absPdg = std::abs(pdg);
    if (!isRelevant(absPdg)) return false;

    // Skip antiparticles unless requested
    if (pdg < 0 && !fReadAnti) return false;

    // Recompute energy from 3-momentum and PDG mass for consistency
    const double mass = massForPdg(pdg);
    const double E    = std::sqrt(px*px + py*py + pz*pz + mass*mass);

    Particle p;
    p.pdg = pdg;
    p.mom.SetPxPyPzE(px, py, pz, E);
    // pos is left at (0,0,0) — will be filled by SourceSampler

    if (absPdg == PDG::kProton) {
        ev.protons.push_back(p);
    } else if (absPdg == PDG::kHe3) {
        ev.he3s.push_back(p);
    }

    return true;
}
