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
#include "SPHParticle.h"
#include "DensityFunctor.h"
#include "HydroForceFunctor.h"
#include "SimpleVtkWriter.h"
#include "SPHConfig.h"

#include "autopas/utils/ArrayMath.h"

using Particle = SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

void Initialize(AutoPasContainer &sphSystem, double density_0) {
  AutoPasLog(INFO, "Initialization started");
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure(density_0);
  }
  AutoPasLog(INFO, "Initialization completed");
}

void velocityVerletFirstStep(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(part->getAcceleration() * dt * 0.5);
    part->addR(part->getV() * dt);
    part->addDensity(part->getDensityDot() * dt);
  }
}

void velocityVerletSecondStep(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(part->getAcceleration() * dt * 0.5);
  }
}

void calculateDensityDot(AutoPasContainer &sphSystem) {
  DensityFunctor<Particle> densityFunctor;

  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setDensityDot(0.);
    densityFunctor.AoSFunctor(*part, *part);
    part->setDensityDot(part->getDensityDot() / 2);
  }

  sphSystem.computeInteractions(&densityFunctor);
}

void updatePressure(AutoPasContainer &sphSystem, double density_0) {
  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure(density_0);
  }
}

void calculateHydroForce(AutoPasContainer &sphSystem, double cutoff, double lj_cutoff,
                         double lj_epsilon, double lj_sigma, double alpha) {
  HydroForceFunctor<Particle> hydroForceFunctor(cutoff, lj_cutoff, lj_epsilon, lj_sigma, alpha);

  AUTOPAS_OPENMP(parallel)
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

void addExternalForce(AutoPasContainer &sphSystem, const std::array<double, 3> &externalForce) {
  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addAcceleration(externalForce);
  }
}

void generateGhostParticles(AutoPasContainer &sphSystem, double cutoff) {
  std::vector<Particle> ghosts;
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();
  bool needs_ghost;
  double min_d;
  double max_d;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    auto pos = part->getR();
    auto vel = part->getV();

    for (size_t dim = 0; dim < 3; dim++) {
      needs_ghost = false;
      min_d = pos[dim] - boxMin[dim];
      max_d = boxMax[dim] - pos[dim];
      if (min_d < cutoff & min_d > 0) {
        needs_ghost = true;
        pos[dim] = -min_d; // Mirrored position
        vel[dim] = -vel[dim]; // Reverse normal velocity (no-slip)
      } else if (max_d < cutoff & max_d > 0) {
        needs_ghost = true;
        pos[dim] = boxMax[dim] + max_d; // Mirrored position
        vel[dim] = -vel[dim]; // Reverse normal velocity (no-slip)
      }
      if (needs_ghost){
        Particle ghost = *part;
        ghost.setR(pos);
        ghost.setV(vel);
        ghost.setIsGhost(true);
        ghosts.push_back(ghost);
      }
    }
  }
  for (auto &g : ghosts) {
    sphSystem.addHaloParticle(g);
  }
}

void addEnteringParticles(AutoPasContainer &sphSystem, std::vector<Particle> &invalidParticles) {
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();

  for (auto &p : invalidParticles) {
    // first we have to correct the position of the particles, s.t. they lie inside of the box.
    auto pos = p.getR();
    auto vel = p.getV();
    for (size_t dim = 0; dim < 3; dim++) {
      if (pos[dim] < boxMin[dim]) {
        // has to be smaller than boxMax
        pos[dim] = std::min(std::nextafter(boxMax[dim], -1), boxMin[dim] + (boxMin[dim] - pos[dim]));
        vel[dim] *= -.9; // -1 would be a perfectly reflective boundary, decimal used as damping
      } else if (pos[dim] >= boxMax[dim]) {
        // should at least be boxMin
        pos[dim] = std::max(boxMin[dim], boxMax[dim] - (pos[dim] - boxMax[dim]));
        vel[dim] *= -.9;
      }
    }
    p.setR(pos);
    p.setV(vel);
    // add moved particles again
    sphSystem.addParticle(p);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " <path_to_config.yaml>\n";
      return 1;
  }

  std::string configFilePath = argv[1];

  SPHConfig config;
  if (!config.loadFromFile(configFilePath)) {
      return 1;
  }

  AutoPasContainer sphSystem;
  std::array<double, 3> boxMin(config.getBoxMin()), boxMax(config.getBoxMax());
  double dt, t_end;
  int write_freq;
  double cutoff, lj_cutoff, lj_epsilon, lj_sigma, density, alpha;
  config.SetupContainer(sphSystem, &dt, &t_end, &write_freq, &cutoff, &lj_cutoff,
                        &lj_epsilon, &lj_sigma, &density, &alpha);

  std::set<autopas::ContainerOption> allowedContainers{autopas::ContainerOption::linkedCells,
                                                       autopas::ContainerOption::verletLists,
                                                       autopas::ContainerOption::verletListsCells};
  sphSystem.setAllowedContainers(allowedContainers);

  std::set<autopas::DataLayoutOption> allowedDataLayouts{autopas::DataLayoutOption::aos};
  sphSystem.setAllowedDataLayouts(allowedDataLayouts);

  sphSystem.init();

  std::array<double, 3> gravity = config.getGravity();
  std::vector<double> forceTimestamps = config.getForceTimestamps();
  std::vector<std::array<double, 3>> customForces = config.getCustomForces();

  std::array<double, 3> externalForce = gravity;

  config.SetupParticles(sphSystem);
  Initialize(sphSystem, density);

  SimpleVtkWriter vtkWriter(config.getSessionName(), config.getOutputFolder(), config.getMaxDigits());

  size_t step = 0;
  size_t force_step = 0;

  for (double time = 0.; time < t_end; time += dt, ++step) {
    velocityVerletFirstStep(sphSystem, dt);

    auto invalidParticles = sphSystem.updateContainer();
    addEnteringParticles(sphSystem, invalidParticles);

    if (time > forceTimestamps[force_step]) {
      externalForce = autopas::utils::ArrayMath::add(gravity, customForces[force_step]);
      force_step += 1;
    }
    config.generateBoundaryParticles(sphSystem);
    generateGhostParticles(sphSystem, cutoff);
    updatePressure(sphSystem, density);
    calculateDensityDot(sphSystem);
    calculateHydroForce(sphSystem, cutoff, lj_cutoff, lj_epsilon, lj_sigma, alpha);
    addExternalForce(sphSystem, externalForce);

    velocityVerletSecondStep(sphSystem, dt);

    if (step % write_freq == 0) {
      AutoPasLog(INFO, "Iteration {} completed", step);
      // AutoPasLog(INFO, "Number of halo particles: {}", sphSystem.getNumberOfParticles(autopas::IteratorBehavior::halo));
      vtkWriter.recordTimestep(step, sphSystem, boxMin, boxMax, autopas::IteratorBehavior::owned);
    }
  }
}

