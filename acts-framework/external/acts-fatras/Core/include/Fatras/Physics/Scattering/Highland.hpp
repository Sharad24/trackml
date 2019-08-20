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

/// @The struct for the Highland-based scattering
///
/// This will scatter particles with a single gaussian distribution
/// according to the highland formula.
struct Highland {

  /// The highland formula
  Acts::detail::HighlandScattering highlandForumla;

  /// @brief Call operator to perform this scattering
  ///
  /// @tparam generator_t is a random number generator type
  /// @tparam detector_t is the detector information type
  /// @tparam particle_t is the particle information type
  ///
  /// @param[in] generator is the random number generator
  /// @param[in] detector the detector information
  /// @param[in] particle the particle which is being scattered
  ///
  /// @return a scattering angle in 3D
  template <typename generator_t, typename detector_t, typename particle_t>
  double operator()(generator_t &generator, const detector_t &detector,
                    particle_t &particle) const {

    // Gauss distribution, will be sampled sampled with generator
    GaussDist gaussDist = GaussDist(0., 1.);

    /// thickness in X0
    double dInX0 = detector.thickness() / detector.material().X0();
    bool electron = std::abs(particle.pdg()) == 11;

    // return projection factor times sigma times grauss radnom
    return M_SQRT2 *
           highlandForumla(particle.p(), particle.beta(), dInX0, electron) *
           gaussDist(generator);
  }
};

} // namespace Fatras
