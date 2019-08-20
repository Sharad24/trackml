// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
class Surface;
}

namespace FW {

namespace Data {

  /// @brief Information to be writtern out per hit surface
  /// - this is the bare information on a surface
  ///
  /// @tparam particle_t
  /// Type of the particle that created to the simhit
  template <typename particle_t>
  struct SimHit
  {
    /// The surface where the hit was created
    const Acts::Surface* surface = nullptr;
    /// The global position of the hit
    Acts::Vector3D position = Acts::Vector3D(0., 0., 0.);
    /// The global direction of the particle at hit position
    Acts::Vector3D direction = Acts::Vector3D(0., 0., 0.);
    /// The value representing the hit (e.g. energy deposit)
    double value = 0.;
    /// The time stamp of the particle
    double timeStamp = 0.;
    /// The particle that created the simulation hit
    particle_t particle;
  };

  /// @brief Constructor of Sensitive Hits, provides a simple
  /// interface to construct simulation hits from Fatras
  struct SimHitCreator
  {

    /// @brief Hit creation from Fatras simulation
    ///
    /// @param surface is the Surface where the hit is created
    /// @param position is the global hit position
    /// @param direction is the momentum direction at hit position
    /// @param value the simulated value
    /// @param time is the timeStamp
    /// @param particle the particle for the truth link
    ///
    /// @return a SimHit object that is created from input parameters
    template <typename particle_t>
    SimHit<particle_t>
    operator()(const Acts::Surface&  surface,
               const Acts::Vector3D& position,
               const Acts::Vector3D& direction,
               double                value,
               double                time,
               const particle_t&     simParticle) const
    {
      SimHit<particle_t> simHit;
      simHit.surface   = &surface;
      simHit.position  = position;
      simHit.direction = direction;
      simHit.value     = value;
      simHit.particle  = simParticle;
      return simHit;
    }
  };

}  // end of namesapce Data
}  // end of namespace FW