// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <mutex>
#include <vector>

#include <Acts/Utilities/Logger.hpp>
#include <Acts/Utilities/Units.hpp>
#include <Pythia8/Pythia.h>

#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"

namespace FW {

class Pythia8Generator
{
public:
  struct Config
  {
    pdg_type pdgBeam0  = 2212;  ///< pdg code of incoming beam 1
    pdg_type pdgBeam1  = 2212;  ///< pdg code of incoming beam 2
    double   cmsEnergy = 14 * Acts::units::_TeV;  ///< center of mass energy
    std::vector<std::string> settings
        = {{"HardQCD:all = on"}};  ///< additional pythia settings
  };

  static std::function<std::vector<Data::SimVertex<Data::SimParticle>>(
      RandomEngine&)>
  makeFunction(const Config& cfg);

  // try to prevent pythia breakage by forbidding copying
  Pythia8Generator()                        = delete;
  Pythia8Generator(const Pythia8Generator&) = delete;
  Pythia8Generator(Pythia8Generator&&)      = default;
  Pythia8Generator&
  operator=(const Pythia8Generator&)
      = delete;
  Pythia8Generator&
  operator=(Pythia8Generator&& other)
      = default;

  Pythia8Generator(const Config&        cfg,
                   Acts::Logging::Level level = Acts::Logging::INFO);

  std::vector<Data::SimVertex<Data::SimParticle>>
  operator()(RandomEngine& rng);

private:
  /// Private access to the logging instance
  const Acts::Logger&
  logger() const
  {
    return (*m_logger);
  }

  Config                              m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  ::Pythia8::Pythia                   m_pythia8;
  std::mutex                          m_pythia8Mutex;
};

}  // namespace FW
