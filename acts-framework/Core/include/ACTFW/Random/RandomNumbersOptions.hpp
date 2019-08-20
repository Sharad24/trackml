// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RandomNumbersOptions.hpp
///////////////////////////////////////////////////////////////////

#ifndef ACTFW_OPTIONS_RANDOMNUMBEROPTIONS_HPP
#define ACTFW_OPTIONS_RANDOMNUMBEROPTIONS_HPP

#include <iostream>
#include "ACTFW/Random/RandomNumbersSvc.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  // common evgen options, with an rnd prefix
  template <class AOPT>
  void
  addRandomNumbersOptions(AOPT& opt)
  {
    opt.add_options()("rnd-seed",
                      po::value<int>()->default_value(1234567890),
                      "Seed of the random number engine.");
  }

  /// read the random number options and return a Config object
  template <class AMAP>
  FW::RandomNumbersSvc::Config
  readRandomNumbersConfig(const AMAP& vm)
  {

    FW::RandomNumbersSvc::Config randomConfig;
    randomConfig.seed = vm["rnd-seed"].template as<int>();
    // return the config
    return randomConfig;
  }
}
}

#endif  // ACTFW_OPTIONS_RANDOMNUMBEROPTIONS_HPP