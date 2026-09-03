#ifndef TOYMCENGINEA4_H
#define TOYMCENGINEA4_H

#include "Config.h"

#include "TDirectory.h"

#include <string>

/**
 * Drievr class for the coalescence toy Monte Carlo engine.
 * It handles the configuration, bookkeeping, and output file management.
 */
class ToyMcEngine {
public:
    ToyMcEngine(const TDirectory* outputFile, std::string& configFile);
    ToyMcEngine(const TDirectory* outputFile, const Config& config) : fOutputFile(const_cast<TDirectory*>(outputFile)), fConfig(config) {};
    void runCoalescence(TDirectory* out = nullptr);
    void runNucleusNucleonSource(TDirectory* out = nullptr);

private:

    TDirectory*     fOutputFile;
    Config          fConfig;

};

#endif // TOYMCENGINEA4_H