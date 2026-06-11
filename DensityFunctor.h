#pragma once

#include "SPHKernels.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"

template <class Particle_T>
class DensityFunctor : public autopas::PairwiseFunctor<Particle_T, DensityFunctor<Particle_T>> {
 public:

  DensityFunctor() : autopas::PairwiseFunctor<Particle_T, DensityFunctor<Particle_T>>(0.){};

  virtual std::string getName() override { return "SPHDensityFunctor"; }

  bool isRelevantForTuning() override { return true; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  /**
   * Calculates the density contribution of the interaction of particle i and j.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  inline void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3 = true) override {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const std::array<double, 3> dr = j.getR() - i.getR();  // ep_j[j].pos - ep_i[i].pos;
    const double density =
        j.getMass() * SPHKernels::W(dr, i.getSmoothingLength());  // ep_j[j].mass * W(dr, ep_i[i].smth)
    i.addDensity(density);
    if (newton3) {
      // Newton 3:
      // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
      const double density2 = i.getMass() * SPHKernels::W(dr, j.getSmoothingLength());
      j.addDensity(density2);
    }
  }
};

