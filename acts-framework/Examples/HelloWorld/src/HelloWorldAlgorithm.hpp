// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_EXAMPLES_HELLOWORLD_H
#define ACTFW_EXAMPLES_HELLOWORLD_H 1

#include <memory>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FWE {

/// @class Algorithm
class HelloWorldAlgorithm : public FW::BareAlgorithm
{
public:
  /// Constructor
  HelloWorldAlgorithm(Acts::Logging::Level level = Acts::Logging::INFO);

  /// Framework execode method
  /// @param [in] context is the Algorithm context for event consistency
  FW::ProcessCode
  execute(FW::AlgorithmContext context) const final override;
};

}  // namespace FWE

#endif  // ACTFW_EXAMPLES_HELLOWORLD_H
