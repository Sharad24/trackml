// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-14
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include <limits>

#include "ACTFW/Framework/BareAlgorithm.hpp"

namespace FW {

/// Select particles by applying some selection cuts.
///
/// Vertices that end up w/o particles are dropped from the output.
class ParticleSelector : public BareAlgorithm
{
public:
  struct Config
  {
    /// The input collection
    std::string input;
    /// The output collection
    std::string output;
    /// Maximum distance from the origin in the transverse plane
    double rhoMax = std::numeric_limits<double>::max();
    /// Maximum absolute distance from the origin along z
    double absZMax = std::numeric_limits<double>::max();
    // Particle cuts
    double phiMin    = std::numeric_limits<double>::lowest();
    double phiMax    = std::numeric_limits<double>::max();
    double etaMin    = std::numeric_limits<double>::lowest();
    double etaMax    = std::numeric_limits<double>::max();
    double absEtaMin = std::numeric_limits<double>::lowest();
    double absEtaMax = std::numeric_limits<double>::max();
    double ptMin     = 0.0;
    double ptMax     = std::numeric_limits<double>::max();
    /// Keep neutral particles
    bool keepNeutral = true;
  };

  ParticleSelector(const Config&        cfg,
                   Acts::Logging::Level level = Acts::Logging::INFO);

  ProcessCode
  execute(AlgorithmContext ctx) const;

private:
  Config m_cfg;
};

}  // namespace FW
