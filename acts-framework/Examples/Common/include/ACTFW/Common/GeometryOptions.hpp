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

  /// the particle gun options, the are prefixes with gp
  template <class AOPT>
  void
  addGeometryOptions(AOPT& opt)
  {
    opt.add_options()("geo-surface-loglevel",
                      po::value<size_t>()->default_value(3),
                      "The outoput log level for the surface building.")(
        "geo-layer-loglevel",
        po::value<size_t>()->default_value(3),
        "The output log level for the layer building.")(
        "geo-volume-loglevel",
        po::value<size_t>()->default_value(3),
        "The output log level for the volume building.")(
        "geo-subdetectors",
        po::value<read_strings>()->multitoken()->default_value({{}}),
        "Sub detectors for the output writing");
  }
}  // namespace Options
}  // namespace FW
