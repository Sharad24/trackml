// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/RootDetector/BuildRootDetector.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "detail/FatrasExampleBase.hpp"

/// @brief adding some specific options for this geometry type
struct RootOptions
{
  /// @brief operator to be called to add options for the generic detector
  ///
  // @tparam options_t Type of the options object
  ///@param opt Options object
  template <typename options_t>
  void
  operator()(options_t& opt)
  {
    FW::Options::addRootGeometryOptions<options_t>(opt);
  }
};

/// @brief geometry getter, the operator() will be called int he example base
struct RootGeometry
{

  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> detElementStore;

  /// @brief operator called to construct the tracking geometry
  ///
  /// @tparam variable_map_t Type of the variable map template for parameters
  ///
  /// @param vm the parameter map object
  ///
  /// @return a closed TrackingGeometry object
  template <typename variable_map_t>
  std::shared_ptr<const Acts::TrackingGeometry>
  operator()(variable_map_t& vm)
  {
    return FW::Root::buildRootDetector<variable_map_t>(vm, detElementStore);
  }
};

/// @brief main executable
///
/// @param argc The argument count
/// @param argv The argument list
int
main(int argc, char* argv[])
{
  // --------------------------------------------------------------------------------
  RootOptions  rootOptions;
  RootGeometry rootGeometry;

  // now process it
  return fatrasExample(argc, argv, rootOptions, rootGeometry);
}