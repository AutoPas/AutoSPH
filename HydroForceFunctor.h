#pragma once

template <class Particle_T>
class HydroForceFunctor : public autopas::PairwiseFunctor<Particle_T, HydroForceFunctor<Particle_T>> {
 private:
  const double _cutoff;
  const double _alpha;
  const double _beta;

 public:

  HydroForceFunctor(double cutoff, double alpha, double beta)
      // the actual cutoff used is dynamic. 0 is used to pass the sanity check.
      : autopas::PairwiseFunctor<Particle_T, HydroForceFunctor<Particle_T>>(cutoff),
        _cutoff{cutoff},
        _alpha{alpha},
        _beta{beta} {};

  virtual std::string getName() override { return "SPHHydroForceFunctor"; }

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

    const double h_ij = 0.5 * (i.getSmoothingLength() + j.getSmoothingLength());
    const double c_ij = 0.5 * (i.getSoundSpeed() + j.getSoundSpeed());
    const double rho_ij = 0.5 * (rho_i + rho_j);
    const double varphi = 0.1 * h_ij;
    const double phi_ij = h_ij * dvdr / (autopas::utils::ArrayMath::dot(dr, dr) + varphi);

    const double AV = (dvdr < 0) ? (-_alpha * c_ij * phi_ij + _beta * phi_ij * phi_ij) / rho_ij : 0;

    const std::array<double, 3> gradW_ij =
        (SPHKernels::gradW(dr, i.getSmoothingLength()) + SPHKernels::gradW(dr, j.getSmoothingLength())) * 0.5;
    // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
    // ep_j[j].smth));

    double scale =
        i.getPressure() / (rho_i * rho_i) + j.getPressure() / (rho_j * rho_j) + AV;
    i.subAcceleration(gradW_ij * (scale * j.getMass()));
    // hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
    // ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) *
    // gradW_ij;
    if (newton3) {
      j.addAcceleration(gradW_ij * (scale * i.getMass()));
      // Newton3, gradW_ij = -gradW_ji
    }

    i.addEngDot(autopas::utils::ArrayMath::dot(gradW_ij, dv) * (scale * j.getMass()));
    // hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
    // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

    if (newton3) {
      j.addEngDot(autopas::utils::ArrayMath::dot(gradW_ij, dv) * (scale * -1 * i.getMass()));
      // Newton 3
    }
  }
};
