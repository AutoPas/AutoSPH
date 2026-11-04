#pragma once

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <string>
#include <vector>

#include "SPHParticle.h"
#include "autopas/AutoPas.h"

using AutoPasContainer = autopas::AutoPas<SPHParticle>;

class SPHConfig {
 private:
  std::string _outputFolder;
  std::string _sessionName;
  unsigned int _maxDigits;
  bool _writeFiles;

  std::array<double, 3> _boxMin{};
  std::array<double, 3> _boxMax{};
  double _timeStep;
  double _totalTime;
  double _writeFrequency;
  double _cutoff;
  double _skinToCutoffRatio;
  unsigned int _rebuildFrequency;
  unsigned int _numSamples;

  std::array<double, 3> _gravity{};
  std::vector<double> _forceTimestamps;
  std::vector<std::array<double, 3>> _customForces;
  std::vector<std::array<double, 3>> _probePositions;

  double _density;
  std::array<double, 3> _particleBoxMin{};
  std::array<double, 3> _particleBoxMax{};
  std::array<unsigned int, 3> _particleNum{};
  double _smoothingLength;
  double _boundarySpacing;
  double _boundaryPressure;
  double _soundSpeed;
  double _alpha;
  double _beta;

  std::array<double, 3> _particleVelocity{0.0, 0.0, 0.0};
  std::array<double, 3> _boundaryParticleSpacing{};
  std::array<double, 3> _particleSpacing;
  double _particlesVolume = 1;
  double _containerVolume = 1;
  double _particleMass;
  unsigned int _totalNumParticles;

 public:
  SPHConfig() = default;

  bool loadFromFile(const std::string &filename) {
    try {
      YAML::Node config = YAML::LoadFile(filename);

      _outputFolder = config["output"]["output_folder"].as<std::string>();
      _sessionName = config["output"]["session_name"].as<std::string>();
      _writeFiles = config["output"]["write_files"].as<bool>(true);

      _boxMin = config["simulation"]["box_min"].as<std::array<double, 3>>();
      _boxMax = config["simulation"]["box_max"].as<std::array<double, 3>>();

      _timeStep = config["simulation"]["time_step"].as<double>();
      _totalTime = config["simulation"]["total_time"].as<double>();
      _writeFrequency = config["simulation"]["write_frequency"].as<double>();

      _skinToCutoffRatio = config["simulation"]["skin_cutoff_ratio"].as<double>();
      _rebuildFrequency = config["simulation"]["rebuild_frequency"].as<unsigned int>();
      _numSamples = config["simulation"]["num_samples"].as<unsigned int>();

      _gravity = config["forces"]["gravity"].as<std::array<double, 3>>();
      if (config["forces"]["time_dependent"]) {
        for (const auto &node : config["forces"]["time_dependent"]) {
          _forceTimestamps.push_back(node["timestamp"].as<double>());
          _customForces.push_back(node["force"].as<std::array<double, 3>>());
        }
      }
      _forceTimestamps.push_back(_totalTime + _timeStep);
      _customForces.push_back({0.0, 0.0, 0.0});

      if (config["probes"]) {
        _probePositions = config["probes"].as<std::vector<std::array<double, 3>>>();
      }

      _density = config["particles"]["density"].as<double>();

      _particleBoxMin = config["particles"]["particle_box_min"].as<std::array<double, 3>>();
      _particleBoxMax = config["particles"]["particle_box_max"].as<std::array<double, 3>>();
      _particleNum = config["particles"]["particle_num"].as<std::array<unsigned int, 3>>();

      _smoothingLength = config["particles"]["smoothing_length"].as<double>();
      _cutoff = SPHKernels::getKernelSupportRadius() * _smoothingLength;

      _boundarySpacing = config["particles"]["boundary_spacing"].as<double>() * _smoothingLength;
      _boundaryPressure = config["particles"]["boundary_pressure"].as<double>();
      _soundSpeed = config["particles"]["sound_speed"].as<double>();
      _alpha = config["particles"]["alpha"].as<double>();
      _beta = config["particles"]["beta"].as<double>();

      return true;
    } catch (const YAML::Exception &e) {
      std::cerr << "Error parsing YAML file: " << e.what() << std::endl;
      return false;
    }
  }

  std::string getOutputFolder() { return _outputFolder; }
  std::string getSessionName() { return _sessionName; }
  unsigned int getMaxDigits() { return std::to_string(static_cast<int>(_totalTime / _timeStep)).length(); }
  std::array<double, 3> getBoxMin() { return _boxMin; }
  std::array<double, 3> getBoxMax() { return _boxMax; }
  std::array<double, 3> getGravity() { return _gravity; }
  std::vector<double> getForceTimestamps() { return _forceTimestamps; }
  std::vector<std::array<double, 3>> getCustomForces() { return _customForces; }
  std::vector<std::array<double, 3>> getProbes() { return _probePositions; }

  void SetupContainer(AutoPasContainer &sphSystem, double *dt, double *t_end, int *writeFreq, bool *writeFiles,
                      double *cutoff, double *smoothingLength, double *density, double *alpha, double *beta,
                      double *boundarySpacing, double *boundaryPressure) {
    sphSystem.setBoxMin(_boxMin);
    sphSystem.setBoxMax(_boxMax);

    sphSystem.setNumSamples(_numSamples);
    sphSystem.setCutoff(_cutoff);
    sphSystem.setVerletSkin(_skinToCutoffRatio * _cutoff);
    sphSystem.setVerletRebuildFrequency(_rebuildFrequency);

    *dt = _timeStep;
    *t_end = _totalTime + _timeStep * 0.5;
    *writeFreq = static_cast<int>(std::round(_writeFrequency / _timeStep));
    *writeFiles = writeFiles;
    *cutoff = _cutoff;
    *smoothingLength = _smoothingLength;
    *density = _density;
    *alpha = _alpha;
    *beta = _beta;
    *boundarySpacing = _boundarySpacing;
    *boundaryPressure = _boundaryPressure;
  }

  void SetupParticles(AutoPasContainer &sphSystem) {
    unsigned int i = 0;
    double num_div;

    AutoPasLog(INFO, "Setup started");

    _totalNumParticles = _particleNum[0] * _particleNum[1] * _particleNum[2];

    for (size_t dim = 0; dim < 3; dim++) {
      num_div = _particleNum[dim] - 1;
      if (_particleBoxMin[dim] <= _boxMin[dim]) {
        _particleBoxMin[dim] = _boxMin[dim];
        num_div += 1;
      }
      if (_particleBoxMax[dim] >= _boxMax[dim]) {
        _particleBoxMax[dim] = _boxMax[dim];
        num_div += 1;
      }

      _particleSpacing[dim] = (_particleBoxMax[dim] - _particleBoxMin[dim]) / num_div;
      _boundaryParticleSpacing[dim] =
          (_boxMax[dim] - _boxMin[dim]) / std::round((_boxMax[dim] - _boxMin[dim]) / _particleSpacing[dim]);

      if (_particleBoxMin[dim] == _boxMin[dim]) {
        _particleBoxMin[dim] += 1 * _particleSpacing[dim];
      }
      if (_particleBoxMax[dim] == _boxMax[dim]) {
        _particleBoxMax[dim] -= 1 * _particleSpacing[dim];
      }

      _particlesVolume *= _particleBoxMax[dim] - _particleBoxMin[dim];
      _containerVolume *= _boxMax[dim] - _boxMin[dim];

      _particleBoxMax[dim] +=
          0.5 *
          _particleSpacing[dim];  // increasing particleBoxMax to ensure particle is added in case of rounding errors
    }

    _particleMass = _particlesVolume * _density / _totalNumParticles;

    for (double x = _particleBoxMin[0]; x < _particleBoxMax[0]; x += _particleSpacing[0]) {
      for (double y = _particleBoxMin[1]; y < _particleBoxMax[1]; y += _particleSpacing[1]) {
        for (double z = _particleBoxMin[2]; z < _particleBoxMax[2]; z += _particleSpacing[2]) {
          SPHParticle ith({x, y, z}, _particleVelocity, i++, _particleMass, _smoothingLength, _soundSpeed);
          ith.setDensity(_density);
          ith.setEnergy(2.5);
          sphSystem.addParticle(ith);
        }
      }
    }

    AutoPasLog(INFO, "Setup completed");

    AutoPasLog(INFO, "Particles Volume: {:.5f}", _particlesVolume);
    AutoPasLog(INFO, "Container Volume: {:.5f}", _containerVolume);
    AutoPasLog(INFO, "Fluid / Container Volume: {:.2f}", _particlesVolume / _containerVolume);

    AutoPasLog(INFO, "Number of particles: {}", sphSystem.getNumberOfParticles());
    AutoPasLog(INFO, "Particle spacing: [{:.5f}, {:.5f}, {:.5f}] ", _particleSpacing[0], _particleSpacing[1],
               _particleSpacing[2]);
    AutoPasLog(INFO, "h / particle spacing: [{:.2f}, {:.2f}, {:.2f}] ", _smoothingLength / _particleSpacing[0],
               _smoothingLength / _particleSpacing[1], _smoothingLength / _particleSpacing[2]);
    AutoPasLog(INFO, "Cutoff distance: {:.3f}", _cutoff);
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
      a_max = _boxMax[dim2] + 0.5 * _boundaryParticleSpacing[dim2];
      b_max = _boxMax[dim3] - 0.5 * _boundaryParticleSpacing[dim3];

      for (double a = _boxMin[dim2]; a < a_max; a += _boundaryParticleSpacing[dim2]) {
        position[dim2] = a;
        for (double b = _boxMin[dim3] + _boundaryParticleSpacing[dim3]; b < b_max;
             b += _boundaryParticleSpacing[dim3]) {
          position[dim3] = b;

          position[dim1] = _boxMin[dim1] - _boundarySpacing;
          SPHParticle p1(position, _particleVelocity, id++, _particleMass, _smoothingLength, _soundSpeed);
          p1.setDensity(_density);
          p1.setIsBoundary(true);
          sphSystem.addHaloParticle(p1);

          position[dim1] = _boxMax[dim1] + _boundarySpacing;
          SPHParticle p2(position, _particleVelocity, id++, _particleMass, _smoothingLength, _soundSpeed);
          p2.setDensity(_density);
          p2.setIsBoundary(true);
          sphSystem.addHaloParticle(p2);
        }
      }
    }
    // 8 corners
    x = std::nextafter(_boxMin[0], _boxMin[0] - 1);
    for (size_t i = 0; i < 2; ++i, x = std::nextafter(_boxMax[0], _boxMax[0] + 1)) {
      y = std::nextafter(_boxMin[1], _boxMin[1] - 1);
      for (size_t j = 0; j < 2; ++j, y = std::nextafter(_boxMax[1], _boxMax[1] + 1)) {
        z = std::nextafter(_boxMin[2], _boxMin[2] - 1);
        for (size_t k = 0; k < 2; ++k, z = std::nextafter(_boxMax[2], _boxMax[2] + 1)) {
          SPHParticle p({x, y, z}, _particleVelocity, id++, _particleMass, _smoothingLength, _soundSpeed);
          p.setDensity(_density);
          p.setIsBoundary(true);
          sphSystem.addHaloParticle(p);
        }
      }
    }
  }
};
