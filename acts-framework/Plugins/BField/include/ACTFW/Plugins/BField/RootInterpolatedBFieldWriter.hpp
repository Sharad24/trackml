// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <TTree.h>
#include <array>
#include <boost/optional.hpp>
#include <mutex>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
class InterpolatedBFieldMap;
}

namespace FW {

namespace BField {

  /// @enum gridType
  /// Describes the axes definition of the grid of the magnetic field map
  enum GridType { rz = 0, xyz = 1 };

  /// @class RootInterpolatedBFieldWriter
  ///
  /// Writes out the Acts::InterpolatedbFieldMap. Currently implemented for 'rz'
  /// and 'xyz' field maps.

  class RootInterpolatedBFieldWriter
  {
  public:
    struct Config
    {
      /// The name of the output tree
      std::string treeName = "TTree";
      /// The name of the output file
      std::string fileName = "TFile.root";
      /// the file access mode (recreate by default)
      std::string fileMode = "recreate";
      /// The magnetic field to be written out
      std::shared_ptr<const Acts::InterpolatedBFieldMap> bField = nullptr;
      /// How the magnetic field map should be written out
      GridType gridType = xyz;
      /// [optional] Setting the range to be printed out in either r (for
      /// cylinder coordinates) or x/y (in cartesian coordinates)
      /// @note setting this parameter is optional, in case no boundaries are
      /// handed over the full magnetic field map will be printed out
      boost::optional<std::array<double, 2>> rBounds;
      /// [optional] Setting the range in z to be printed out
      /// @note setting this parameter is optional, in case no boundaries are
      /// handed over the full magnetic field map will be printed out
      boost::optional<std::array<double, 2>> zBounds;
      /// Number of bins in r
      /// @note setting this parameter is optional, in case no bin numbers are
      /// handed over the full magnetic field map will be printed out
      size_t rBins = 200;
      /// Number of bins in z
      // @note setting this parameter is optional, in case no bin numbers are
      /// handed over the full magnetic field map will be printed out
      size_t zBins = 300;
      /// Number of bins in phi
      // @note setting this parameter is optional, in case no bin numbers are
      /// handed over the full magnetic field map will be printed out
      size_t phiBins = 100;
    };

    /// Write down an interpolated magnetic field map
    static void
    run(const Config&                       cfg,
        std::unique_ptr<const Acts::Logger> logger
        = Acts::getDefaultLogger("RootInterpolatedBFieldWriter",
                                 Acts::Logging::INFO));
  };

}  // namespace BField

}  // namespace FW
