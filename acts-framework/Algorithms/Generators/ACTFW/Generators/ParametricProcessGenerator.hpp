// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cmath>
#include <vector>

#include <Acts/Utilities/Units.hpp>

#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"

namespace FW {

/// Generate particles from uniform parameter distributions.
///
/// Generates a single process vertex with the given number of tracks. Each
/// track's momentum and direction is drawn from uniform parameter
/// distributions. Position and time are always set to zero.
class ParametricProcessGenerator
{
public:
  struct Config
  {
    /// Number of particles
    size_t numParticles = 1;
    /// Low, high for the transverse point of closest approach
    std::array<double, 2> d0Range = {{0.0, 0.0}};
    /// Low, high for the z position at the point of closest approach
    std::array<double, 2> z0Range = {{0.0, 0.0}};
    /// Low, high for the transverse angle
    std::array<double, 2> phiRange = {{-M_PI, M_PI}};
    /// Low, high for pseudo-rapidity
    std::array<double, 2> etaRange = {{-4.0, 4.0}};
    /// Low, high for transverse momentum
    std::array<double, 2> ptRange
        = {{100 * Acts::units::_MeV, 10 * Acts::units::_GeV}};
    /// Particle mass
    double mass = 0.;
    /// Particle charge
    double charge = 1.0;
    /// Randomize the charge; this will also flip the pdg identifier sign
    bool randomCharge = false;
    /// (Absolute) pdg type of the particle
    pdg_type pdg = 0;
  };

  ParametricProcessGenerator(const Config& cfg);

  /// Generate a single process vertex with the given number of particles.
  std::vector<Data::SimVertex<Data::SimParticle>>
  operator()(RandomEngine& rng) const;

private:
  Config m_cfg;
};

}  // namespace FW
