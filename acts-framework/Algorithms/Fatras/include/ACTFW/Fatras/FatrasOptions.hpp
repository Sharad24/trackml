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
#include "Fatras/Kernel/SelectorList.hpp"
#include "Fatras/Selectors/ChargeSelectors.hpp"
#include "Fatras/Selectors/KinematicCasts.hpp"
#include "Fatras/Selectors/SelectorHelpers.hpp"
#include "FatrasAlgorithm.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  /// @brief read the Fatras options
  ///
  /// Adding Fatras specific options to the Options package
  ///
  /// @tparam aopt_t Type of the options object (API bound to boost)
  ///
  /// @param [in] opt_t The options object where the specific digitization
  /// options are attached to
  template <typename aopt_t>
  void
  addFatrasOptions(aopt_t& opt)
  {
    opt.add_options()(
        "fatras-sim-particles",
        po::value<std::string>()->default_value("fatras-particles"),
        "The collection of simulated particles.")(
        "fatras-sim-hits",
        po::value<std::string>()->default_value("fatras-hits"),
        "The collection of simulated hits")(
        "fatras-em-ionisation",
        po::value<bool>()->default_value(true),
        "Switch on ionisiation loss of charged particles")(
        "fatras-em-radiation",
        po::value<bool>()->default_value(true),
        "Switch on radiation for charged particles")(
        "fatras-em-scattering",
        po::value<bool>()->default_value(true),
        "Switch on multiple scattering")(
        "fatras-em-conversions",
        po::value<bool>()->default_value(false),
        "Switch on gamma conversions")("fatras-had-interaction",
                                       po::value<bool>()->default_value(false),
                                       "Switch on hadronic interaction")(
        "fatras-debug-output",
        po::value<bool>()->default_value(false),
        "Switch on debug output on/off");
  }

  /// @brief read the fatras specific options and return a Config file
  ///
  ///@tparam omap_t Type of the options map
  ///@param vm the options map to be read out
  template <typename AMAP,
            typename simulator_t,
            typename event_collection_t,
            typename hit_t>
  typename FatrasAlgorithm<simulator_t, event_collection_t, hit_t>::Config
  readFatrasConfig(const AMAP& vm, simulator_t& simulator)
  {
    // Create a config
    typename FatrasAlgorithm<simulator_t, event_collection_t, hit_t>::Config
        fatrasConfig(std::move(simulator));

    // set the collections
    fatrasConfig.simulatedHitCollection
        = vm["fatras-sim-hits"].template as<std::string>();
    fatrasConfig.simulatedEventCollection
        = vm["fatras-sim-particles"].template as<std::string>();

    typedef Fatras::ChargedSelector            CSelector;
    typedef Fatras::Max<Fatras::casts::absEta> CMaxEtaAbs;
    typedef Fatras::Min<Fatras::casts::pT>     CMinPt;
    typedef Fatras::SelectorListAND<CSelector, CMinPt, CMaxEtaAbs>
        ChargedSelector;

    typedef Fatras::NeutralSelector            NSelector;
    typedef Fatras::Max<Fatras::casts::absEta> NMaxEtaAbs;
    typedef Fatras::Min<Fatras::casts::E>      NMinE;

    simulator.chargedSelector.template get<CMaxEtaAbs>().valMax = 5.;
    simulator.chargedSelector.template get<CMinPt>().valMin
        = 100. * Acts::units::_MeV;

    simulator.neutralSelector.template get<NMaxEtaAbs>().valMax = 5.;
    simulator.neutralSelector.template get<NMinE>().valMin
        = 100. * Acts::units::_MeV;

    // and return the config
    return fatrasConfig;
  }

}  // namespace Options
}  // namespace FW
