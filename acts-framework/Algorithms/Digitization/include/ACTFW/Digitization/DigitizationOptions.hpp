// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Units.hpp"
#include "DigitizationAlgorithm.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  /// @brief Digitization options
  /// Adds specific digitization options to the boost::program_options
  ///
  /// @tparam aopt_t Type of the options object (API bound to boost)
  ///
  /// @param [in] opt_t The options object where the specific digitization
  /// options are attached to
  template <typename aopt_t>
  void
  addDigitizationOptions(aopt_t& opt)
  {
    opt.add_options()("digi-spacepoints",
                      po::value<std::string>()->default_value("space-points"),
                      "Collection name of the produced space points.")(
        "digi-clusters",
        po::value<std::string>()->default_value("clusters"),
        "Collection name of the produced clustes.")(
        "digi-resolution-file",
        po::value<std::string>()->default_value(""),
        "Name of the resolution file (root format).");
  }

  ///@brief  Read the digitization options and return a Config object
  ///
  ///@tparam omap_t Type of the options map
  ///@param vm the options map to be read out
  template <typename omap_t>
  DigitizationAlgorithm::Config
  readDigitizationConfig(const omap_t& vm)
  {
    // create a config
    DigitizationAlgorithm::Config digiConfig;
    digiConfig.spacePointCollection
        = vm["digi-spacepoints"].template as<std::string>();
    digiConfig.clusterCollection
        = vm["digi-clusters"].template as<std::string>();
    digiConfig.resolutionFile
        = vm["digi-resolution-file"].template as<std::string>();
    // and return the config
    return digiConfig;
  }
}  // namespace Options
}  // namespace FW
