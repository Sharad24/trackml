// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <vector>
#include "ACTFW/Barcode/Barcode.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace FW {

// Typedef the process code
typedef unsigned int process_code;

namespace Data {

  /// @brief Vertex information struct for physics process samplers:
  /// - all quatities are calculated at first construction as they may
  ///   be used by downstream samplers
  ///
  /// @note if a sampler changes one of the parameters, consistency
  /// can be broken, so it should update the rest (no checking done)
  template <typename particle_t = SimParticle>
  struct SimVertex
  {

    /// The vertex position
    Acts::Vector3D position = Acts::Vector3D(0., 0., 0.);

    /// The ingoing particles in the vertex
    std::vector<particle_t> in = {};

    /// The outgoing particles from the vertex
    std::vector<particle_t> out = {};

    /// An optional process code
    process_code processCode = 9;

    /// An optional time stamp
    double timeStamp = 0.;

    /// Default
    SimVertex() = default;

    /// @brief Construct a particle consistently
    ///
    /// @param ertex The vertex position
    /// @param in The ingoing particles - copy
    /// @param out The outgoing particles (copy - can we do a move ?)
    /// @param vprocess The process code
    /// @param time The time stamp of this vertex
    SimVertex(const Acts::Vector3D&          vertex,
              const std::vector<particle_t>& ingoing  = {},
              std::vector<particle_t>        outgoing = {},
              process_code                   process  = 0,
              double                         time     = 0.)
      : position(vertex)
      , in(ingoing)
      , out(outgoing)
      , processCode(process)
      , timeStamp(time)
    {
    }

    /// Forward the particle access to the outgoing particles: begin
    ///
    /// @tparam particle_t Type of the particle
    typename std::vector<particle_t>::iterator
    outgoing_begin()
    {
      return out.begin();
    }

    /// Forward the particle access to the outgoing particles: end
    ///
    /// @tparam particle_t Type of the particle
    typename std::vector<particle_t>::iterator
    outgoing_end()
    {
      return out.end();
    }

    // Outgoing particles
    const std::vector<particle_t>&
    outgoing() const
    {
      return out;
    }

    /// Forward the particle access to the outgoing particles: insert
    ///
    /// @tparam particle_t Type of the particle
    ///
    /// @param inparticles are the particles to be inserted
    typename std::vector<particle_t>::iterator
    outgoing_insert(const std::vector<particle_t>& inparticles)
    {
      return out.insert(out.end(), inparticles.begin(), inparticles.end());
    }
  };

}  // end of namespace Data
}  // end of namespace FW