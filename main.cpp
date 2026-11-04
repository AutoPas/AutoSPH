#include <mpi.h>
#include <array>
#include <iostream>

#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"

#include "SPHParticle.h"
#include "DensityForceFunctor.h"
#include "SimpleVtkWriter.h"
#include "SPHConfig.h"
#include "TerminalOutput.h"

using Particle = SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

void InitializePressure(AutoPasContainer &sphSystem, double density_0) {
  AutoPasLog(INFO, "Pressure initialization started");
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure(density_0);
  }
  AutoPasLog(INFO, "Pressure initialization completed");
}

void eulerStep(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(part->getAcceleration() * dt);
    part->addR(part->getV() * dt);
    part->addDensity(part->getDensityDot() * dt);
    part->addEnergy(part->getEngDot() * dt);
  }
}

void velocityVerletFirstStep(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(part->getAcceleration() * dt * 0.5);
    part->addR(part->getV() * dt);
    part->addDensity(part->getDensityDot() * dt);
    part->addEnergy(part->getEngDot() * dt);
  }
}

void velocityVerletSecondStep(AutoPasContainer &sphSystem, const double dt) {
  using namespace autopas::utils::ArrayMath::literals;

  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addV(part->getAcceleration() * dt * 0.5);
  }
}

void updatePressure(AutoPasContainer &sphSystem, double density_0) {
  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure(density_0);
  }
}

void calculateDensityForce(AutoPasContainer &sphSystem, DensityForceFunctor<Particle> &densityForceFunctor) {
  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->setDensityDot(0.);
    part->setAcceleration(std::array<double, 3>{0., 0., 0.});
    part->setEngDot(0.);
  }
  sphSystem.computeInteractions(&densityForceFunctor);
}

void addExternalForce(AutoPasContainer &sphSystem, const std::array<double, 3> &externalForce) {
  AUTOPAS_OPENMP(parallel)
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addAcceleration(externalForce);
  }
}

bool calculateGhostPosVel(double &pos, double &vel, double boxMin, double boxMax, double cutoff, double boundarySpacing) {
  // If the particle is close enough to a boundary, it returns true and
  // calculates the position and velocity of the ghost particle in 1 dimension.
  double min_d = pos - boxMin;
  double max_d = boxMax - pos;

  if (min_d < cutoff) {
    pos = boxMin - min_d - boundarySpacing;
    vel = -vel;
    return true;
  } else if (max_d < cutoff) {
    pos = boxMax + max_d + boundarySpacing;
    vel = -vel;
    return true;
  }
  return false;
}

bool addGhostPosVel(size_t dim, double cutoff, double boundarySpacing, std::array<double, 3> &boxMin, std::array<double, 3> &boxMax,
                    std::vector<std::array<double, 3>> &positions, std::vector<std::array<double, 3>> &velocities, 
                    std::array<double, 3> &newPos, std::array<double, 3> &newVel) {
  double pos = positions[0][dim];
  double vel = velocities[0][dim];

  if (calculateGhostPosVel(pos, vel, boxMin[dim], boxMax[dim], cutoff, boundarySpacing)){
    newPos[dim] = pos;
    newVel[dim] = vel;
    positions.push_back(newPos);
    velocities.push_back(newVel);
    return true;
  }
  return false;
}

void generateGhostParticles(AutoPasContainer &sphSystem, double cutoff, double boundarySpacing, double boundaryPressure) {
  std::array<double, 3> boxMin = sphSystem.getBoxMin();
  std::array<double, 3> boxMax = sphSystem.getBoxMax();
  std::array<double, 3> newPos, newVel;
  std::vector<Particle> ghosts;

  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    std::vector<std::array<double, 3>> positions = { part->getR() };
    std::vector<std::array<double, 3>> velocities = { part->getV() };

    // Checks if the particle is close to a boundary in each dimension and creates 1 ghost if it is within cutoff distance
    // of 1 wall; creates 3 ghosts in case it is close to 2 walls; and 7 ghosts if it is close to 3 walls.
    for (size_t dim1 = 0; dim1 < 3; dim1++) {
      newPos = positions[0];
      newVel = velocities[0];
      if (addGhostPosVel(dim1, cutoff, boundarySpacing, boxMin, boxMax, positions, velocities, newPos, newVel)) {
        for (size_t dim2 = dim1 + 1; dim2 < 3; dim2++) {
          if (addGhostPosVel(dim2, cutoff, boundarySpacing, boxMin, boxMax, positions, velocities, newPos, newVel)) {
            if (dim2 == 1) {
              size_t dim3 = 2;
              addGhostPosVel(dim3, cutoff, boundarySpacing, boxMin, boxMax, positions, velocities, newPos, newVel);
              newPos = positions[1];
              newVel = velocities[1];
            }
          }
        }
      }
    }
    for (size_t i = 1; i < positions.size(); ++i) {
        Particle ghost = *part;
        ghost.setR(positions[i]);
        ghost.setV(velocities[i]);
        ghost.setPressure(part->getPressure() * boundaryPressure);
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
    auto pos = p.getR();
    auto vel = p.getV();
    for (size_t dim = 0; dim < 3; dim++) {
      if (pos[dim] < boxMin[dim]) {
        pos[dim] = std::min(std::nextafter(boxMax[dim], -1), boxMin[dim] + (boxMin[dim] - pos[dim]));
        vel[dim] *= -1;
      } else if (pos[dim] >= boxMax[dim]) {
        pos[dim] = std::max(boxMin[dim], boxMax[dim] - (pos[dim] - boxMax[dim]));
        vel[dim] *= -1;
      }
    }
    p.setR(pos);
    p.setV(vel);
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
  double dt, t_end;
  int writeFreq;
  bool writeFiles;
  double cutoff, smoothingLength, density, alpha, beta, boundarySpacing, boundaryPressure;
  config.SetupContainer(sphSystem, &dt, &t_end, &writeFreq, &writeFiles, &cutoff, &smoothingLength,
                        &density, &alpha, &beta, &boundarySpacing, &boundaryPressure);

  std::set<autopas::ContainerOption> allowedContainers{autopas::ContainerOption::directSum,
                                                       autopas::ContainerOption::linkedCells,
                                                       autopas::ContainerOption::linkedCellsReferences,
                                                       autopas::ContainerOption::verletLists,
                                                       autopas::ContainerOption::verletListsCells,
                                                       autopas::ContainerOption::verletClusterLists,
                                                       autopas::ContainerOption::varVerletListsAsBuild,
                                                       autopas::ContainerOption::pairwiseVerletLists,
                                                       autopas::ContainerOption::octree};
  sphSystem.setAllowedContainers(allowedContainers);

  std::set<autopas::DataLayoutOption> allowedDataLayouts{autopas::DataLayoutOption::aos};
  sphSystem.setAllowedDataLayouts(allowedDataLayouts);

  sphSystem.init();

  std::array<double, 3> boxMin(config.getBoxMin()), boxMax(config.getBoxMax());
  std::array<double, 3> gravity = config.getGravity();
  std::array<double, 3> externalForce = gravity;
  std::vector<double> forceTimestamps = config.getForceTimestamps();
  std::vector<std::array<double, 3>> customForces = config.getCustomForces();
  std::vector<std::array<double, 3>> probes = config.getProbes();

  config.SetupParticles(sphSystem);
  InitializePressure(sphSystem, density);

  SimpleVtkWriter vtkWriter(config.getSessionName(), config.getOutputFolder(), config.getMaxDigits(), probes.size(), configFilePath);
  TerminalOutput terminalOutput;

  DensityForceFunctor<Particle> densityForceFunctor(cutoff, alpha, beta);

  size_t step = 0;
  size_t force_step = 0;
  size_t fileStep = 0;
  const size_t maxIterations = static_cast<size_t>(t_end/dt);

  AutoPasLog(INFO, "Simulation started");

  for (double time = 0.; time < t_end; time += dt, ++step) {
    // eulerStep(sphSystem, dt);
    velocityVerletFirstStep(sphSystem, dt);

    auto invalidParticles = sphSystem.updateContainer();
    addEnteringParticles(sphSystem, invalidParticles);

    if (time > forceTimestamps[force_step]) {
      externalForce = autopas::utils::ArrayMath::add(gravity, customForces[force_step]);
      force_step += 1;
    }
    config.generateBoundaryParticles(sphSystem);
    updatePressure(sphSystem, density);
    generateGhostParticles(sphSystem, cutoff, boundarySpacing, boundaryPressure);
    calculateDensityForce(sphSystem, densityForceFunctor);
    addExternalForce(sphSystem, externalForce);

    velocityVerletSecondStep(sphSystem, dt);

    if (step % writeFreq == 0) {
      terminalOutput.printProgress(step, maxIterations);
      if (writeFiles){
        vtkWriter.recordTimestep(step, sphSystem, boxMin, boxMax, autopas::IteratorBehavior::owned);
        vtkWriter.recordProbes(step, time, fileStep, cutoff, smoothingLength, sphSystem, probes);
        ++fileStep;
      }
    }
  }

  AutoPasLog(INFO, "Simulation completed");
}

