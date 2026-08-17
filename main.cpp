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
#include "SPHParticle.cpp"
#include "DensityFunctor.h"
#include "HydroForceFunctor.h"
#include "SimpleVtkWriter.h"

using Particle = SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

void SetupIC(AutoPasContainer &sphSystem, double *dt, double *end_time, double density, const std::array<double, 3> &bBoxMax) {
  // Place SPH particles
  AutoPasLog(INFO, "Setup started");

  const double i_box = 0.5;
  const int part_num = 16;
  const double part_cube_size = bBoxMax[0] * i_box;
  const double dx = bBoxMax[2] / (part_num + 1);
  const double dim0 = bBoxMax[0] - 2 * dx;
  const double dim1 = bBoxMax[1] * 0.5 - dx - std::fmod(bBoxMax[0] * 0.5, dx);
  const double dim2 = bBoxMax[2] - 2 * dx;
  const double total_mass =  dim0 * dim1 * dim2 * density;
  const double part_mass = total_mass / (part_num * part_num * part_num * 0.5);
  unsigned int i = 0;
  for (double x = dx; x < bBoxMax[0]; x += dx) {         // NOLINT
    // for (double y = bBoxMax[1] - bBoxMax[1] * i_box; y < bBoxMax[1] - dx; y += dx) {       // NOLINT
    for (double y = dx; y < bBoxMax[1]*0.5; y += dx) {       // NOLINT
      for (double z = dx; z < bBoxMax[2]; z += dx) {     // NOLINT
        Particle ith({x, y, z}, {0, 0, 0}, i++, part_mass, 0.012, 20.0);
        ith.setDensity(density);
        ith.setEnergy(2.5);
        sphSystem.addParticle(ith);
      }
    }
  }

  // Set dt and end time
  *dt = .0002;
  *end_time = 5;

  AutoPasLog(INFO, "Setup completed");
  AutoPasLog(INFO, "Number of particles (i): {}", i);
  AutoPasLog(INFO, "Number of particles: {}", sphSystem.getNumberOfParticles());
}

void Initialize(AutoPasContainer &sphSystem, double density_0) {
  AutoPasLog(INFO, "Initialization started");
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure(density_0);
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

void updatePressure(AutoPasContainer &sphSystem, double density_0) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->calcPressure(density_0);
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

void addGravity(AutoPasContainer &sphSystem, const std::array<double, 3> &gravity) {
  for (auto part = sphSystem.begin(autopas::IteratorBehavior::owned); part.isValid(); ++part) {
    part->addAcceleration(gravity);
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

    for (auto dim = 0; dim < 3; dim++) {
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
        ghost.setIsBoundary(true);
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
    for (auto dim = 0; dim < 3; dim++) {
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

int main() {
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = boxMax[1] = boxMax[2] = .25;
  double cutoff = 0.03;               // 0.012*2.5=0.03; where 2.5 = kernel support radius
  unsigned int rebuildFrequency = 6;  // has to be multiple of two, as there are two functor calls per iteration.
  double skinToCutoffRatio = 0.15;
  std::array<double, 3> gravity({0., -10., 0});
  double slosh_acc = 5;
  double density = 1000.0;

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
  SetupIC(sphSystem, &dt, &t_end, density, boxMax);
  Initialize(sphSystem, density);
  const int record_freq = static_cast<int>(std::round(0.005 / dt));
  // LogParticlePositions(sphSystem);

  SimpleVtkWriter vtkWriter("serial_test_run", "./output", 5);

  applyConstantForce(sphSystem);
  size_t step = 0;
  for (double time = 0.; time < t_end; time += dt, ++step) {
    generateGhostParticles(sphSystem, cutoff);

    calculateDensity(sphSystem);
    updatePressure(sphSystem, density);
    calculateHydroForce(sphSystem);
    addGravity(sphSystem, gravity);

    eulerStep(sphSystem, dt);

    if (step % record_freq == 0) {
      AutoPasLog(INFO, "Iteration {} completed", step);
      AutoPasLog(INFO, "Number of particles: {}", sphSystem.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo));
      vtkWriter.recordTimestep(step, sphSystem, boxMin, boxMax);

      if (time > t_end * .8) {
        gravity = {0., -10., 0};
      } else if (time > t_end * .6) {
        gravity = {slosh_acc, -10., 0};
      } else if (time > t_end * .4) {
        gravity = {0., -10., 0};
      } else if (time > t_end * .2) {
        gravity = {-slosh_acc, -10., 0};
      }
    }

    auto invalidParticles = sphSystem.updateContainer();
    addEnteringParticles(sphSystem, invalidParticles);
  }

  // LogParticlePositions(sphSystem);
}

