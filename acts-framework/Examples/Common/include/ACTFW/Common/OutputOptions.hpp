// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdlib>
#include <iostream>
#include <utility>
#include "ACTFW/Utilities/Options.hpp"

namespace po = boost::program_options;

namespace FW {
namespace Options {

  /// @brief the common output options that are added to the
  /// job options
  ///
  /// @tparam aopt_t Type of the options object, bound to boost API
  /// @param [in] opt The options object for attaching specific options
  template <typename aopt_t>
  void
  addOutputOptions(aopt_t& opt)
  {
    // Add specific options for this example
    opt.add_options()("output-dir",
                      po::value<std::string>()->default_value(""),
                      "Output directory location.")(
        "output-root",
        po::value<bool>()->default_value(false),
        "Switch on to write '.root' output file(s).")(
        "output-csv",
        po::value<bool>()->default_value(false),
        "Switch on to write '.csv' output file(s).")(
        "output-obj",
        po::value<bool>()->default_value(false),
        "Switch on to write '.obj' ouput file(s).")(
        "output-json",
        po::value<bool>()->default_value(false),
        "Switch on to write '.json' ouput file(s).");
  }
}  // namespace Options
}  // namespace FW