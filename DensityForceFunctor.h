#pragma once

#include "SPHKernels.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"

template <class Particle_T>
class DensityForceFunctor : public autopas::PairwiseFunctor<Particle_T, DensityForceFunctor<Particle_T>> {
 private:
  const double _cutoff;
  const double _alpha;
  const double _beta;

 public:

  DensityForceFunctor(double cutoff, double alpha, double beta)
      // the actual cutoff used is dynamic. 0 is used to pass the sanity check.
      : autopas::PairwiseFunctor<Particle_T, DensityForceFunctor<Particle_T>>(cutoff),
        _cutoff{cutoff},
        _alpha{alpha},
        _beta{beta} {};

  virtual std::string getName() override { return "SPHDensityForceFunctor"; }

  bool isRelevantForTuning() override { return true; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

  /**
   * Calculates the contribution of the interaction of particle i and j to the
   * hydrodynamic force.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3 = true) override {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }

    const std::array<double, 3> dr = i.getR() - j.getR();
    // const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;

    const double distance = autopas::utils::ArrayMath::L2Norm(dr);

    if (distance >= _cutoff) {
      return;
    }

    const std::array<double, 3> dv = i.getV() - j.getV();
    // const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

    double dvdr = autopas::utils::ArrayMath::dot(dv, dr);

    const double rho_i = i.getDensity();
    const double rho_j = j.getDensity();
    const double mass_i = i.getMass();
    const double mass_j = j.getMass();

    const double h_ij = 0.5 * (i.getSmoothingLength() + j.getSmoothingLength());
    const double c_ij = 0.5 * (i.getSoundSpeed() + j.getSoundSpeed());
    const double rho_ij = 0.5 * (rho_i + rho_j);
    const double varphi = 0.1 * h_ij;
    const double phi_ij = h_ij * dvdr / (autopas::utils::ArrayMath::dot(dr, dr) + varphi);

    const double AV = (dvdr < 0) ? (-_alpha * c_ij * phi_ij + _beta * phi_ij * phi_ij) / rho_ij : 0;

    const std::array<double, 3> gradW_ij = SPHKernels::gradW(dr, h_ij);
    const double dv_gradW = autopas::utils::ArrayMath::dot(dv, gradW_ij);
    const double scale = i.getPressure() / (rho_i * rho_i) + j.getPressure() / (rho_j * rho_j) + AV;

    i.addDensityDot(mass_j * dv_gradW);
    i.subAcceleration(gradW_ij * (scale * mass_j));
    i.addEngDot(autopas::utils::ArrayMath::dot(gradW_ij, dv) * (scale * mass_j));

    if (newton3) {
      j.addDensityDot(mass_i * dv_gradW);
      j.addAcceleration(gradW_ij * (scale * mass_i));
      j.addEngDot(autopas::utils::ArrayMath::dot(gradW_ij, dv) * (scale * -1 * mass_i));
    }
  }
};
