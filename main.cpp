/**
 * @file sph-main-mpi.cpp
 * @date 10.04.2018
 * @author seckler
 */

#include <mpi.h>

#include <array>
#include <cmath>
#include <iostream>

#include "autopas/AutoPas.h"
// #include "autopas/particles/ParticleDefinitions.h"
#include "SPHParticle.cpp"

using Particle = SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

// using AutoPasContainer = autopas::AutoPas<autopas::ParticleBaseFP64>;

void SetupIC(AutoPasContainer &sphSystem, double *end_time, const std::array<double, 3> &bBoxMax) {
  // Place SPH particles
  AutoPasLog(INFO, "Setup started");

  const double dx = 1.0 / 2.0;
  unsigned int i = 0;
  for (double x = 0; x < bBoxMax[0]; x += dx) {         // NOLINT
    for (double y = 0; y < bBoxMax[1]; y += dx) {       // NOLINT
      for (double z = 0; z < bBoxMax[2]; z += dx) {     // NOLINT
        Particle ith({x, y, z}, {0, 0, 0}, i++, 0.75, 0.012, 0.);
        ith.setDensity(1.0);
        ith.setEnergy(2.5);
        sphSystem.addParticle(ith);
      }
    }
  }

  // Set the end time
  *end_time = .012;

  AutoPasLog(INFO, "Setup completed");
  AutoPasLog(INFO, "Number of particles: {}", i);
}

void Initialize(AutoPasContainer &sphSystem) {
  AutoPasLog(INFO, "Initialization started");
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure();
  }
  AutoPasLog(INFO, "Initialization completed");
}

void LogParticlePositions(AutoPasContainer &sphSystem) {
  std::cout << "initialize... completed" << std::endl;
  std::array<double, 3> position;
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    // std::cout << part->getMass() << std::endl;
    position = part->getR();
    AutoPasLog(INFO, "Position of particle {}: {}, {}, {}", part->getID(), position[0], position[1], position[2]);
  }
  std::cout << "initialize... completed" << std::endl;
}

int main() {
  std::cout<< "Hello from SPH world!\n";

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = boxMax[1] = boxMax[2] = 1.;
  double cutoff = 0.03;               // 0.012*2.5=0.03; where 2.5 = kernel support radius
  unsigned int rebuildFrequency = 6;  // has to be multiple of two, as there are two functor calls per iteration.
  double skinToCutoffRatio = 0.15;

  AutoPasContainer sphSystem;
  sphSystem.setNumSamples(
      6);  // has to be multiple of 2, should also be multiple of rebuildFrequency (but this is not necessary).
  sphSystem.setBoxMin(boxMin);
  sphSystem.setBoxMax(boxMax);
  sphSystem.setCutoff(cutoff);
  sphSystem.setVerletSkin(skinToCutoffRatio * cutoff);
  sphSystem.setVerletRebuildFrequency(rebuildFrequency);

  std::set<autopas::ContainerOption> allowedContainers{autopas::ContainerOption::linkedCells,
                                                       autopas::ContainerOption::verletLists,
                                                       autopas::ContainerOption::verletListsCells};
  sphSystem.setAllowedContainers(allowedContainers);

  sphSystem.init();

  double dt;
  double t_end;
  SetupIC(sphSystem, &t_end, boxMax);
  Initialize(sphSystem);
  LogParticlePositions(sphSystem);
}

