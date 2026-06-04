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

int main() {

  std::cout<< "Hello from SPH world!\n";

  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 1.;
  boxMax[1] = boxMax[2] = boxMax[0] / 8.0;
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
}

