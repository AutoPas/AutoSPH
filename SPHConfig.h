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
    double timeStep;
    double totalTime;
    double boxMaxX;
    double boxMaxY;
    double boxMaxZ;
    double cutoff;
    double skinToCutoffRatio;
    unsigned int rebuildFrequency;
    unsigned int numSamples;

public:
    SPHConfig() = default;

    bool loadFromFile(const std::string& filename) {
        try {
            YAML::Node config = YAML::LoadFile(filename);

            timeStep = config["simulation"]["time_step"].as<double>();
            totalTime = config["simulation"]["total_time"].as<double>();

            boxMaxX = config["simulation"]["BoxSize_X"].as<double>();
            boxMaxY = config["simulation"]["BoxSize_Y"].as<double>();
            boxMaxZ = config["simulation"]["BoxSize_Z"].as<double>();

            cutoff = config["simulation"]["cutoff"].as<double>();
            skinToCutoffRatio = config["simulation"]["skin_cutoff_ratio"].as<double>();
            rebuildFrequency = config["simulation"]["rebuild_frequency"].as<unsigned int>();
            numSamples = config["simulation"]["num_samples"].as<unsigned int>();

            return true;
        } catch (const YAML::Exception& e) {
            std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
            return false;
        }
    }

    double getTimeStep() const { return timeStep; }
    double getTotalTime() const { return totalTime; }
    double getBoxMaxX() const { return boxMaxX; }
    double getBoxMaxY() const { return boxMaxY; }
    double getBoxMaxZ() const { return boxMaxZ; }

    void SetupContainer(AutoPasContainer &sphSystem) {
        sphSystem.setNumSamples(numSamples);
        sphSystem.setCutoff(cutoff);
        sphSystem.setVerletSkin(skinToCutoffRatio * cutoff);
        sphSystem.setVerletRebuildFrequency(rebuildFrequency);
    }
};
