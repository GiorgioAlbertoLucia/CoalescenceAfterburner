#include "ToyMcEngine.h"
#include "CoalescenceEngineFactory.h"
#include "BookKeeping.h"

#include <iostream>
#include <cmath>

ToyMcEngine::ToyMcEngine(const TDirectory* outputFile, std::string& configFile):
    fOutputFile(const_cast<TDirectory*>(outputFile))
{
    fConfig.loadFromFile(configFile);
}

void ToyMcEngine::run(TDirectory* out) {
    
    if (!out) {
        out = fOutputFile;
    }
    BookKeeping bookKeeping;
    
    // create the appropriate coalescence engine based on the nucleus name
    CoalescenceEngine* engine = CoalescenceEngineFactory::createCoalescenceEngine(fConfig);
    engine->run(bookKeeping);
    
    bookKeeping.write(out);
}
