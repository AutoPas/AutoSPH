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
    std::array<double, 3> boundaryParticleSpacing{};
    double particleMass;
    unsigned int total_num_particles;

    double lj_cutoff;
    double lj_epsilon;
    double lj_sigma;

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

            lj_cutoff = config["LJ potential"]["lj_cutoff"].as<double>();
            lj_epsilon = config["LJ potential"]["lj_epsilon"].as<double>();
            lj_sigma = config["LJ potential"]["lj_sigma"].as<double>();

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

    void SetupContainer(AutoPasContainer &sphSystem, double *dt, double *t_end, int *write_freq,
                        double *_cutoff, double *_lj_cutoff, double *_lj_epsilon, double *_lj_sigma,
                        double *_density, double *_alpha) {
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
        *_lj_cutoff = lj_cutoff;
        *_lj_epsilon = lj_epsilon;
        *_lj_sigma = lj_sigma;
        *_density = density;
        *_alpha = alpha;
    }

    void SetupParticles(AutoPasContainer &sphSystem) {
        unsigned int i = 0;
        double num_div;
        double particlesVolume = 1;
        std::array<double, 3> particleSpacing;

        AutoPasLog(INFO, "Setup started");

        total_num_particles = particleNum[0] * particleNum[1] * particleNum[2];

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
            boundaryParticleSpacing[dim] = (boxMax[dim] - boxMin[dim]) / std::round((boxMax[dim] - boxMin[dim]) / particleSpacing[dim]);

            if (particleBoxMin[dim] == boxMin[dim]) { particleBoxMin[dim] += 1 * particleSpacing[dim]; }
            if (particleBoxMax[dim] == boxMax[dim]) { particleBoxMax[dim] -= 1 * particleSpacing[dim]; }

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

    void generateBoundaryParticles(AutoPasContainer &sphSystem) {
        std::array<double, 3> position{};
        size_t dim1, dim2, dim3;
        int id = 0;
        double a_max, b_max;
        double x, y, z;

        for (dim1 = 0; dim1 < 3; dim1++) {
            dim2 = (dim1 + 1) % 3;
            dim3 = (dim1 + 2) % 3;
            a_max = boxMax[dim2] + 0.5 * boundaryParticleSpacing[dim2];
            b_max = boxMax[dim3] - 0.5 * boundaryParticleSpacing[dim3];

            for (double a = boxMin[dim2]; a < a_max; a += boundaryParticleSpacing[dim2]) {
                position[dim2] = a;
                for (double b = boxMin[dim3] + boundaryParticleSpacing[dim3]; b < b_max;
                     b += boundaryParticleSpacing[dim3]) {
                    position[dim3] = b;
                    position[dim1] = std::nextafter(boxMin[dim1], boxMin[dim1] - 1);
                    for (size_t i = 0; i < 2; ++i, position[dim1] = std::nextafter(boxMax[dim1], boxMax[dim1] + 1)) {
                        SPHParticle p(position, particleVelocity, id++, particleMass, smoothingLength, soundSpeed);
                        p.setDensity(density);
                        p.setIsBoundary(true);
                        sphSystem.addHaloParticle(p);
                    }
                }
            }
        }
        // 8 corners
        x = std::nextafter(boxMin[0], boxMin[0] - 1);
        for (size_t i = 0; i < 2; ++i, x = std::nextafter(boxMax[0], boxMax[0] + 1)) {
            y = std::nextafter(boxMin[1], boxMin[1] - 1);
            for (size_t j = 0; j < 2; ++j, y = std::nextafter(boxMax[1], boxMax[1] + 1)) {
                z = std::nextafter(boxMin[2], boxMin[2] - 1);
                for (size_t k = 0; k < 2; ++k, z = std::nextafter(boxMax[2], boxMax[2] + 1)) {
                    SPHParticle p({x, y, z}, particleVelocity, id++, particleMass, smoothingLength, soundSpeed);
                    p.setDensity(density);
                    p.setIsBoundary(true);
                    sphSystem.addHaloParticle(p);
                }
            }
        }
    }
};
