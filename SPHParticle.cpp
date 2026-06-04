
#pragma once

#include <cstring>
#include <vector>

#include "autopas/particles/ParticleDefinitions.h"

/**
 * Basic SPHParticle class.
 */
class SPHParticle : public autopas::ParticleBaseFP64 {
 public:
  /**
   * Default constructor of SPHParticle.
   * Will initialize all values to some basic defaults.
   */
  SPHParticle()
      : autopas::ParticleBaseFP64(),
        _density(0.),
        _pressure(0.),
        _mass(0.),
        _smth(0.),
        _snds(0.),
        // temporaries / helpers
        _v_sig_max(0.),
        _acc{0., 0., 0.},
        _energy_dot(0.),
        _energy(0.),
        _dt(0.),
        _vel_half{0., 0., 0.},
        _eng_half(0.) {}
 private:
  double _density;   // density
  double _pressure;  // pressure
  double _mass;      // mass
  double _smth;      // smoothing length
  double _snds;      // speed of sound

  // temporaries / helpers
  double _v_sig_max;

  // integrator
  std::array<double, 3> _acc;  // acceleration
  double _energy_dot;          // time derivative of the energy
  double _energy;              // energy
  double _dt;                  // local timestep allowed by this particle

  // integrator
  std::array<double, 3> _vel_half;  // velocity at half time-step
  double _eng_half;                 // energy at half time-step
};
