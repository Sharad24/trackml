// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"

#include "Fatras/Kernel/detail/RandomNumberDistributions.hpp"

namespace Fatras {

/// @brief The struct for the EnergyLoss physics list
///
/// This generates the energy loss according to the Bethe-Bloch
/// description, applying a landau generated enery loss
///
/// It follows the interface of EnergyLoss samplers in Fatras
/// that could return radiated photons for further processing,
/// however, for the Bethe-Bloch application the return vector
/// is always 0.
struct BetheBloch {

  /// Scaling for most probable value
  double scaleFactorMPV = 1.;

  /// Scaling for Sigma
  double scaleFactorSigma = 1.;

  /// Bethe-Bloch Ionisation loss formula
  Acts::detail::IonisationLoss ionisationLoss;

  /// @brief Call operator for the Bethe Bloch energy loss
  ///
  /// @tparam generator_t is a random number generator type
  /// @tparam detector_t is the detector information type
  /// @tparam particle_t is the particle information type
  ///
  /// @param[in] generator is the random number generator
  /// @param[in] detector the detector information
  /// @param[in] particle the particle which is being scattered
  ///
  /// @return empty vector for BetheBloch - no secondaries created
  template <typename generator_t, typename detector_t, typename particle_t>
  std::vector<particle_t> operator()(generator_t &generator,
                                     const detector_t &detector,
                                     particle_t &particle) const {

    // Create a random landau distribution between in the intervall [0,1]
    LandauDist landauDist = LandauDist(0., 1.);
    double landau = landauDist(generator);

    auto eLoss =
        ionisationLoss.dEds(particle.m(), particle.beta(), particle.gamma(),
                            detector.material(), detector.thickness(), false);
    // the actual energy loss
    double energyLoss = eLoss.first;
    // the uncertainty of the mean energy loss
    double energyLossSigma = eLoss.second;

    // Simulate the energy loss
    double deltaE = scaleFactorMPV * std::fabs(energyLoss) +
                    scaleFactorSigma * energyLossSigma * landau;
    // apply the energy loss
    particle.energyLoss(energyLoss);

    // return empty children
    return {};
  }
};

} // namespace Fatras
