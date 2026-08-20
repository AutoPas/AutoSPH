#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "autopas/AutoPas.h"
#include "SPHParticle.h"

using AutoPasContainer = autopas::AutoPas<SPHParticle>;

class SPHConfig {
private:
    std::array<double, 3> boxMin;
    std::array<double, 3> boxMax;
    double timeStep;
    double totalTime;
    double writeFrequency;
    double cutoff;
    double skinToCutoffRatio;
    unsigned int rebuildFrequency;
    unsigned int numSamples;
    std::array<double, 3> gravity;
    std::vector<double> forceTimestamps;
    std::vector<std::array<double, 3>> customForces;

public:
    SPHConfig() = default;

    bool loadFromFile(const std::string& filename) {
        try {
            YAML::Node config = YAML::LoadFile(filename);

            boxMin = config["simulation"]["box_min"].as<std::array<double, 3>>();
            boxMax = config["simulation"]["box_max"].as<std::array<double, 3>>();

            timeStep = config["simulation"]["time_step"].as<double>();
            totalTime = config["simulation"]["total_time"].as<double>();
            writeFrequency = config["simulation"]["write_frequency"].as<double>();

            cutoff = config["simulation"]["cutoff"].as<double>();
            skinToCutoffRatio = config["simulation"]["skin_cutoff_ratio"].as<double>();
            rebuildFrequency = config["simulation"]["rebuild_frequency"].as<unsigned int>();
            numSamples = config["simulation"]["num_samples"].as<unsigned int>();

            gravity = config["forces"]["gravity"].as<std::array<double, 3>>();
            if (config["forces"]["time_dependent"]) {
                for (const auto& node : config["forces"]["time_dependent"]) {
                    forceTimestamps.push_back(node["timestamp"].as<double>());
                    customForces.push_back(node["force"].as<std::array<double, 3>>());
                }
            }
            return true;
        } catch (const YAML::Exception& e) {
            std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
            return false;
        }
    }

    std::array<double, 3> getBoxMin() { return boxMin; }
    std::array<double, 3> getBoxMax() { return boxMax; }
    std::array<double, 3> getGravity() { return gravity; }
    std::vector<double> getForceTimestamps() { return forceTimestamps; }
    std::vector<std::array<double, 3>> getCustomForces() { return customForces; }

    void SetupContainer(AutoPasContainer &sphSystem, double *dt, double *t_end, int *write_freq, double *_cutoff) {
        sphSystem.setBoxMin(boxMin);
        sphSystem.setBoxMax(boxMax);

        sphSystem.setNumSamples(numSamples);
        sphSystem.setCutoff(cutoff);
        sphSystem.setVerletSkin(skinToCutoffRatio * cutoff);
        sphSystem.setVerletRebuildFrequency(rebuildFrequency);

        *dt = timeStep;
        *t_end = totalTime;
        *write_freq = static_cast<int>(std::round(writeFrequency / timeStep));
        *_cutoff = cutoff;
    }
};
