// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_EXAMPLES_RANDOMNUMBERSALGORITHM_H
#define ACTFW_EXAMPLES_RANDOMNUMBERSALGORITHM_H 1

#include <array>
#include <memory>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {
class WhiteBoard;
class RandomNumbersSvc;
}  // namespace FW

namespace FWE {

class RandomNumbersAlgorithm : public FW::BareAlgorithm
{
public:
  struct Config
  {
    std::shared_ptr<FW::RandomNumbersSvc> randomNumbers = nullptr;

    std::array<double, 2> gaussParameters   = {{0., 1.}};
    std::array<double, 2> uniformParameters = {{0., 1.}};
    std::array<double, 2> landauParameters  = {{0., 1.}};
    std::array<double, 2> gammaParameters   = {{0., 1.}};
    int poissonParameter = 40;

    size_t drawsPerEvent = 0;
  };

  /// Constructor
  RandomNumbersAlgorithm(const Config&        cnf,
                         Acts::Logging::Level level = Acts::Logging::INFO);

  /// Framework execode method
  FW::ProcessCode
  execute(FW::AlgorithmContext context) const final override;

private:
  Config m_cfg;
};

}  // namespace FWE

#endif  // ACTFW_EXAMPLES_RANDOMNUMBERSALGORITHM_H
