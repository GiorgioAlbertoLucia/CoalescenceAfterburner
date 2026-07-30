#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

struct Config {

    std::string nucleusName;
    double sourceRadius;
    std::string inputPtHistogramFile;
    std::string inputPtHistogramName;
    std::string outputFile;

    int nEvents = static_cast<int>(1e7);
    int nThreads = 1;
    int randomSeed = 42;

    void loadFromFile(const std::string& filename) {
        YAML::Node config = YAML::LoadFile(filename);

        if (!config["nucleusName"] || !config["sourceRadius"] || !config["inputPtHistogramFile"] || !config["inputPtHistogramName"] || !config["outputFile"]) {
            throw std::invalid_argument("Config: missing required configuration parameters in " + filename);
        }

        nucleusName = config["nucleusName"].as<std::string>();
        sourceRadius = config["sourceRadius"].as<double>();
        inputPtHistogramFile = config["inputPtHistogramFile"].as<std::string>();
        inputPtHistogramName = config["inputPtHistogramName"].as<std::string>();
        outputFile = config["outputFile"].as<std::string>();
        if (config["nEvents"]) {
            nEvents = config["nEvents"].as<int>();
        }
        if (config["nThreads"]) {
            nThreads = config["nThreads"].as<int>();
        }
        if (config["randomSeed"]) {
            randomSeed = config["randomSeed"].as<int>();
        }
    }

};

#endif // CONFIG_H