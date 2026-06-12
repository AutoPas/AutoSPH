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
#include "DensityFunctor.h"
#include "HydroForceFunctor.h"

using Particle = SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

// using AutoPasContainer = autopas::AutoPas<autopas::ParticleBaseFP64>;

void SetupIC(AutoPasContainer &sphSystem, double *dt, double *end_time, const std::array<double, 3> &bBoxMax) {
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

  // Set dt and end time
  *dt = .002;
  *end_time = .024;

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
  std::array<double, 3> position;
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    position = part->getR();
    AutoPasLog(INFO, "Position of particle {}: {}, {}, {}", part->getID(), position[0], position[1], position[2]);
  }
}

void eulerStep(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(part->getAcceleration() * dt);
    part->addR(part->getV() * dt);
  }
}

void applyConstantForce(AutoPasContainer &sphSystem) {
  using namespace autopas::utils::ArrayMath::literals;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setAcceleration({2e3, 1e3, 0.0});
  }
}

void calculateDensity(AutoPasContainer &sphSystem) {
  DensityFunctor<Particle> densityFunctor;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setDensity(0.);
    densityFunctor.AoSFunctor(*part, *part);
    part->setDensity(part->getDensity() / 2);
  }

  sphSystem.computeInteractions(&densityFunctor);
}

void updatePressure(AutoPasContainer &sphSystem) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure();
  }
}

void calculateHydroForce(AutoPasContainer &sphSystem) {
  HydroForceFunctor<Particle> hydroForceFunctor;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    // self interaction leeds to:
    // 1) vsigmax = 2*part->getSoundSpeed()
    // 2) no change in acceleration
    part->setVSigMax(2 * part->getSoundSpeed());
    part->setAcceleration(std::array<double, 3>{0., 0., 0.});
    part->setEngDot(0.);
  }

  sphSystem.computeInteractions(&hydroForceFunctor);
}

void addEnteringParticles(AutoPasContainer &sphSystem, std::vector<Particle> &invalidParticles) {
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();

  for (auto &p : invalidParticles) {
    // first we have to correct the position of the particles, s.t. they lie inside of the box.
    auto pos = p.getR();
    for (auto dim = 0; dim < 3; dim++) {
      if (pos[dim] < boxMin[dim]) {
        // has to be smaller than boxMax
        pos[dim] = std::min(std::nextafter(boxMax[dim], -1), boxMin[dim] + (boxMin[dim] - pos[dim]));
      } else if (pos[dim] >= boxMax[dim]) {
        // should at least be boxMin
        pos[dim] = std::max(boxMin[dim], boxMax[dim] - (pos[dim] - boxMax[dim]));
      }
    }
    p.setR(pos);
    // add moved particles again
    sphSystem.addParticle(p);
  }
}

int main() {
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

  std::set<autopas::DataLayoutOption> allowedDataLayouts{autopas::DataLayoutOption::aos};
  sphSystem.setAllowedDataLayouts(allowedDataLayouts);

  sphSystem.init();

  double dt;
  double t_end;
  SetupIC(sphSystem, &dt, &t_end, boxMax);
  Initialize(sphSystem);
  LogParticlePositions(sphSystem);

  applyConstantForce(sphSystem);
  size_t step = 0;
  for (double time = 0.; time < t_end; time += dt, ++step) {
    calculateDensity(sphSystem);
    updatePressure(sphSystem);
    calculateHydroForce(sphSystem);

    eulerStep(sphSystem, dt);

    auto invalidParticles = sphSystem.updateContainer();
    addEnteringParticles(sphSystem, invalidParticles);

    AutoPasLog(INFO, "Iteration {} completed", step);
  }

  LogParticlePositions(sphSystem);
}

