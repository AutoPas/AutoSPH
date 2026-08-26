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
    std::string outputFolder;
    std::string sessionName;
    unsigned int maxDigits;

    std::array<double, 3> boxMin{};
    std::array<double, 3> boxMax{};
    double timeStep;
    double totalTime;
    double writeFrequency;
    double cutoff;
    double skinToCutoffRatio;
    unsigned int rebuildFrequency;
    unsigned int numSamples;

    std::array<double, 3> gravity{};
    std::vector<double> forceTimestamps;
    std::vector<std::array<double, 3>> customForces;

    double density;
    std::array<double, 3> particleBoxMin{};
    std::array<double, 3> particleBoxMax{};
    std::array<unsigned int, 3> particleNum{};
    double smoothingLength;
    double soundSpeed;
    double alpha;
    std::array<double, 3> particleVelocity{0.0, 0.0, 0.0};

public:
    SPHConfig() = default;

    bool loadFromFile(const std::string& filename) {
        try {
            YAML::Node config = YAML::LoadFile(filename);

            outputFolder = config["output"]["output_folder"].as<std::string>();
            sessionName = config["output"]["session_name"].as<std::string>();

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
            forceTimestamps.push_back(totalTime + timeStep);
            customForces.push_back({0.0, 0.0, 0.0});

            density = config["particles"]["density"].as<double>();

            particleBoxMin = config["particles"]["particle_box_min"].as<std::array<double, 3>>();
            particleBoxMax = config["particles"]["particle_box_max"].as<std::array<double, 3>>();
            particleNum = config["particles"]["particle_num"].as<std::array<unsigned int, 3>>();

            smoothingLength = config["particles"]["smoothing_length"].as<double>();
            soundSpeed = config["particles"]["sound_speed"].as<double>();
            alpha = config["particles"]["alpha"].as<double>();

            return true;
        } catch (const YAML::Exception& e) {
            std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
            return false;
        }
    }

    std::string getOutputFolder() { return outputFolder; }
    std::string getSessionName() { return sessionName; }
    unsigned int getMaxDigits() { return std::to_string(static_cast<int>(totalTime / timeStep)).length(); }
    std::array<double, 3> getBoxMin() { return boxMin; }
    std::array<double, 3> getBoxMax() { return boxMax; }
    std::array<double, 3> getGravity() { return gravity; }
    std::vector<double> getForceTimestamps() { return forceTimestamps; }
    std::vector<std::array<double, 3>> getCustomForces() { return customForces; }

    void SetupContainer(AutoPasContainer &sphSystem, double *dt, double *t_end, int *write_freq, double *_cutoff, double *_density, double *_alpha) {
        sphSystem.setBoxMin(boxMin);
        sphSystem.setBoxMax(boxMax);

        sphSystem.setNumSamples(numSamples);
        sphSystem.setCutoff(cutoff);
        sphSystem.setVerletSkin(skinToCutoffRatio * cutoff);
        sphSystem.setVerletRebuildFrequency(rebuildFrequency);

        *dt = timeStep;
        *t_end = totalTime + timeStep * 0.5;
        *write_freq = static_cast<int>(std::round(writeFrequency / timeStep));
        *_cutoff = cutoff;
        *_density = density;
        *_alpha = alpha;
    }

    void SetupParticles(AutoPasContainer &sphSystem) {
        unsigned int num_div;
        unsigned int i = 0;
        unsigned int total_num_particles = particleNum[0] * particleNum[1] * particleNum[2];
        double x, y, z;
        double particlesVolume = 1;
        double particleMass;
        std::array<double, 3> particleSpacing;

        AutoPasLog(INFO, "Setup started");

        for (size_t dim = 0; dim < 3; dim++) {
            num_div = particleNum[dim] - 1;
            if (particleBoxMin[dim] <= boxMin[dim]) {
                particleBoxMin[dim] = boxMin[dim];
                num_div += 1;
            }
            if (particleBoxMax[dim] >= boxMax[dim]) {
                particleBoxMax[dim] = boxMax[dim];
                num_div += 1;
            }

            particleSpacing[dim] = (particleBoxMax[dim] - particleBoxMin[dim]) / num_div;

            if (particleBoxMin[dim] == boxMin[dim]) { particleBoxMin[dim] += particleSpacing[dim]; }
            if (particleBoxMax[dim] == boxMax[dim]) { particleBoxMax[dim] -= particleSpacing[dim]; }

            particlesVolume *= particleBoxMax[dim] - particleBoxMin[dim];

            particleBoxMax[dim] += 0.5 * particleSpacing[dim]; // increasing particleBoxMax to ensure particle is added in case of rounding errors
        }

        particleMass = particlesVolume * density / total_num_particles;

        for (double x = particleBoxMin[0]; x < particleBoxMax[0]; x += particleSpacing[0]) {
            for (double y = particleBoxMin[1]; y < particleBoxMax[1]; y += particleSpacing[1]) {
                for (double z = particleBoxMin[2]; z < particleBoxMax[2]; z += particleSpacing[2]) {
                    SPHParticle ith({x, y, z}, particleVelocity, i++, particleMass, smoothingLength, soundSpeed);
                    ith.setDensity(density);
                    ith.setEnergy(2.5);
                    sphSystem.addParticle(ith);
                }
            }
        }

        AutoPasLog(INFO, "Setup completed");
        AutoPasLog(INFO, "Number of particles (i): {}", i);
        AutoPasLog(INFO, "Number of particles: {}", sphSystem.getNumberOfParticles());
    }
};
