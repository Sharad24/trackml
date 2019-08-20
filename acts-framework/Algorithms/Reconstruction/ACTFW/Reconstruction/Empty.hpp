// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include "ACTFW/Framework/BareAlgorithm.hpp"

namespace FW {

/// @brief An empty reconstruction algorithm to be use as a template.
class EmptyReconstructionAlgorithm : public BareAlgorithm
{
public:
  struct Config
  {
    /// Input space point collection
    std::string spacePointCollection;
    /// Some config option that changes the reconstruction
    double aConfigOption = 1.23;
  };

  /// Construct the algorithm.
  ///
  /// @param [in] cfg is the configuration struct
  /// @param [in] loglevel is the logging level
  EmptyReconstructionAlgorithm(const Config&        cfg,
                               Acts::Logging::Level logLevel);

  /// Execute the algorithm and reconstruct nothing.
  ///
  /// @param [in] ctx is the algorithm context for event consistency
  /// @return is a process code indicating success or not
  FW::ProcessCode
  execute(AlgorithmContext ctx) const final override;

private:
  Config m_cfg;
};

}  // namespace FW
