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
#include "TerminalOutput.h"

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

bool calculateGhostPosVel(double &pos, double &vel, double boxMin, double boxMax, double cutoff) {
  double min_d = pos - boxMin;
  double max_d = boxMax - pos;

  if (min_d < cutoff & min_d > 0) {
    pos = -min_d; // Mirrored position
    vel = -vel; // Reverse normal velocity (no-slip)
    return true;
  } else if (max_d < cutoff & max_d > 0) {
    pos = boxMax + max_d; // Mirrored position
    vel = -vel; // Reverse normal velocity (no-slip)
    return true;
  }
  return false;
}


void generateGhostParticles(AutoPasContainer &sphSystem, double cutoff) {
  std::vector<Particle> ghosts;
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();
  std::array<double, 3> newPosition, newVelocity;
  bool needs_ghost;
  double pos, vel;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    std::vector<std::array<double, 3>> positions = { part->getR() };
    std::vector<std::array<double, 3>> velocities = { part->getV() };
    for (size_t dim = 0; dim < 3; dim++) {
      pos = positions[0][dim];
      vel = velocities[0][dim];
      if (calculateGhostPosVel(pos, vel, boxMin[dim], boxMax[dim], cutoff)){
        newPosition = positions[0];
        newVelocity = velocities[0];
        newPosition[dim] = pos;
        newVelocity[dim] = vel;
        positions.push_back(newPosition);
        velocities.push_back(newVelocity);
        for (size_t dim2 = dim + 1; dim2 < 3; dim2++) {
          pos = positions[0][dim2];
          vel = velocities[0][dim2];
          if (calculateGhostPosVel(pos, vel, boxMin[dim2], boxMax[dim2], cutoff)){
            newPosition[dim2] = pos;
            newVelocity[dim2] = vel;
            positions.push_back(newPosition);
            velocities.push_back(newVelocity);
            for (size_t dim3 = dim2 + 1; dim3 < 3; dim3++) {
              pos = positions[0][dim3];
              vel = velocities[0][dim3];
              if (calculateGhostPosVel(pos, vel, boxMin[dim3], boxMax[dim3], cutoff)){
                newPosition[dim3] = pos;
                newVelocity[dim3] = vel;
                positions.push_back(newPosition);
                velocities.push_back(newVelocity);
                newPosition = positions[1];
                newVelocity = velocities[1];
              }
            }
          }
        }
      }
    }
    for (size_t i = 1; i < positions.size(); ++i) {
        Particle ghost = *part;
        ghost.setR(positions[i]);
        ghost.setV(velocities[i]);
        ghost.setIsGhost(true);
        ghosts.push_back(ghost);
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
        vel[dim] *= -1;
      } else if (pos[dim] >= boxMax[dim]) {
        // should at least be boxMin
        pos[dim] = std::max(boxMin[dim], boxMax[dim] - (pos[dim] - boxMax[dim]));
        vel[dim] *= -1;
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
  TerminalOutput terminalOutput;

  size_t step = 0;
  size_t force_step = 0;
  const size_t maxIterations = static_cast<size_t>(t_end/dt);

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
      // AutoPasLog(INFO, "Iteration {} completed", step);
      // AutoPasLog(INFO, "Number of halo particles: {}", sphSystem.getNumberOfParticles(autopas::IteratorBehavior::halo));
      terminalOutput.printProgress(step, maxIterations);
      vtkWriter.recordTimestep(step, sphSystem, boxMin, boxMax, autopas::IteratorBehavior::owned);
    }
  }

  AutoPasLog(INFO, "Simulation completed");
}

