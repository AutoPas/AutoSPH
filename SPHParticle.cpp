
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

  /**
   * Constructor of the SPHParticle class.
   * @param r position of the particle
   * @param v velocity of the particle
   * @param id id of the particle. This id should be unique
   * @param mass mass of the particle
   * @param smth smoothing length of the particle
   * @param snds speed of sound (SouND Speed)
   */
  SPHParticle(const std::array<double, 3> &r, const std::array<double, 3> &v, unsigned long id, double mass,
              double smth, double snds)
      : autopas::ParticleBaseFP64(r, v, id),
        _density(0.),
        _pressure(0.),
        _mass(mass),
        _smth(smth),
        _snds(snds),
        // temporaries / helpers
        _v_sig_max(0.),
        _acc{0., 0., 0.},
        _energy_dot(0.),
        _energy(0.),
        _dt(0.),
        _vel_half{0., 0., 0.},
        _eng_half(0.) {}

  /**
   * Destructor of the SPHParticle
   */
  ~SPHParticle() override = default;

  /**
   * Attribute names for the soa arrays
   */
  enum AttributeNames : int {
    ptr,
    mass,
    posX,
    posY,
    posZ,
    smth,
    density,
    velX,
    velY,
    velZ,
    soundSpeed,
    pressure,
    vsigmax,
    accX,
    accY,
    accZ,
    engDot,
    ownershipState
  };

  double getMass() const { return _mass; }

  double getSmoothingLength() const { return _smth; }

  double getDensity() const { return _density; }

  double getSoundSpeed() const { return _snds; }

  double getPressure() const { return _pressure; }

  double getVSigMax() const { return _v_sig_max; }

  const std::array<double, 3> &getAcceleration() const { return _acc; }

  const std::array<double, 3> &getVel_half() const { return _vel_half; }

  double getEngDot() const { return _energy_dot; }

  void setAcceleration(const std::array<double, 3> &acc) { _acc = acc; }

  void setDensity(double density) { _density = density; }

  void setEnergy(double energy) { _energy = energy; }

  void setVel_half(const std::array<double, 3> &vel_half) { SPHParticle::_vel_half = vel_half; }

  void calcPressure() {
    const double hcr = 1.4;
    _pressure = (hcr - 1.0) * _density * _energy;
    _snds = sqrt(hcr * _pressure / _density);
  }

  /**
   * SoA arrays type, cf. AttributeNames
   */
  using SoAArraysType = autopas::utils::SoAType<SPHParticle *,
                                                double,  // mass
                                                double,  // posX
                                                double,  // posY
                                                double,  // posZ
                                                double,  // smth
                                                double,  // density
                                                double,  // velX
                                                double,  // velY
                                                double,  // velZ
                                                double,  // soundSpeed
                                                double,  // pressure
                                                double,  // vsigmax
                                                double,  // accX
                                                double,  // accY
                                                double,  // accZ
                                                double,  // engDot
                                                autopas::OwnershipState>::Type;

  /**
   * Non-const getter for the pointer of this object.
   * @tparam attribute Attribute name.
   * @return this.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   */
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == AttributeNames::mass) {
      return getMass();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::smth) {
      return getSmoothingLength();
    } else if constexpr (attribute == AttributeNames::density) {
      return getDensity();
    } else if constexpr (attribute == AttributeNames::velX) {
      return getV()[0];
    } else if constexpr (attribute == AttributeNames::velY) {
      return getV()[1];
    } else if constexpr (attribute == AttributeNames::velZ) {
      return getV()[2];
    } else if constexpr (attribute == AttributeNames::soundSpeed) {
      return getSoundSpeed();
    } else if constexpr (attribute == AttributeNames::pressure) {
      return getPressure();
    } else if constexpr (attribute == AttributeNames::vsigmax) {
      return getVSigMax();
    } else if constexpr (attribute == AttributeNames::accX) {
      return getAcceleration()[0];
    } else if constexpr (attribute == AttributeNames::accY) {
      return getAcceleration()[1];
    } else if constexpr (attribute == AttributeNames::accZ) {
      return getAcceleration()[2];
    } else if constexpr (attribute == AttributeNames::engDot) {
      return getEngDot();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("SPHParticle::get: unknown attribute");
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   */
  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == AttributeNames::mass) {
      setMass(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      _r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      _r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      _r[2] = value;
    } else if constexpr (attribute == AttributeNames::smth) {
      setSmoothingLength(value);
    } else if constexpr (attribute == AttributeNames::density) {
      setDensity(value);
    } else if constexpr (attribute == AttributeNames::velX) {
      _v[0] = value;
    } else if constexpr (attribute == AttributeNames::velY) {
      _v[1] = value;
    } else if constexpr (attribute == AttributeNames::velZ) {
      _v[2] = value;
    } else if constexpr (attribute == AttributeNames::soundSpeed) {
      setSoundSpeed(value);
    } else if constexpr (attribute == AttributeNames::pressure) {
      setPressure(value);
    } else if constexpr (attribute == AttributeNames::vsigmax) {
      setVSigMax(value);
    } else if constexpr (attribute == AttributeNames::accX) {
      _acc[0] = value;
    } else if constexpr (attribute == AttributeNames::accY) {
      _acc[1] = value;
    } else if constexpr (attribute == AttributeNames::accZ) {
      _acc[2] = value;
    } else if constexpr (attribute == AttributeNames::engDot) {
      setEngDot(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      _ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("SPHParticle::set: unknown attribute");
    }
  }

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
