#ifndef COALESCENCEENGINEFACTORY_H
#define COALESCENCEENGINEFACTORY_H

#include "CoalescenceEngine.h"
#include "CoalescenceEngineHe3.h"
#include "CoalescenceEngineHe4.h"

class CoalescenceEngineFactory {
public:
    static CoalescenceEngine* createCoalescenceEngine(Config& config);
};

CoalescenceEngine* CoalescenceEngineFactory::createCoalescenceEngine(Config& config) {
    if (config.nucleusName == "He3") {
        return new CoalescenceEngineHe3(config);
    } else if (config.nucleusName == "He4") {
        return new CoalescenceEngineHe4(config);
    } else {
        throw std::runtime_error("Unsupported nucleus type: " + config.nucleusName);
    }
}

#endif