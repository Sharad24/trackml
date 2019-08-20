// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CommonOptions.hpp
///////////////////////////////////////////////////////////////////

#pragma once

#include <utility>
#include "Acts/Utilities/Logger.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  /// Add common options
  ///
  /// @tparam Type of the boost option object
  ///
  /// @param defaultEvent is the number of event to be executed
  /// @param defaultValue is the log level default value
  template <typename aopt_t>
  void
  addCommonOptions(aopt_t& opt)
  {
    opt.add_options()("help", "Produce help message")(
        "events,n",
        po::value<size_t>()->default_value(1),
        "The number of events to be processed")(
        "loglevel,l",
        po::value<size_t>()->default_value(2),
        "The output log level. Please set the wished number (0 = VERBOSE, 1 = "
        "DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).");
  }

  /// Read standard options: log level to be returned
  ///
  /// @tparam amap_t Type of the options map
  ///
  /// @param[in] vm Map to be read in
  template <typename amap_t>
  Acts::Logging::Level
  readLogLevel(const amap_t& vm)
  {
    Acts::Logging::Level logLevel
        = Acts::Logging::Level(vm["loglevel"].template as<size_t>());
    logLevel = Acts::Logging::Level(vm["loglevel"].template as<size_t>());
    return logLevel;
  }

  /// Read standard options: number of events to be processed
  ///
  /// @tparam amap_t Type of the options map
  ///
  /// @param[in] vm Map to be read in
  template <typename amap_t>
  size_t
  readNumberOfEvents(const amap_t& vm)
  {
    size_t nEvents = vm["events"].template as<size_t>();
    return nEvents;
  }
}
}
