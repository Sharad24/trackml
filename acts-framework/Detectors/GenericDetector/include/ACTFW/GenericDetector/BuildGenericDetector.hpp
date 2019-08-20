// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
class TrackingGeometry;
}

namespace FW {

namespace Generic {

  /// Global method to build the generic tracking geometry
  /// @param lvl is the debug logging level
  /// @param version is the detector version
  /// return a unique vector to the tracking geometry
  std::unique_ptr<const Acts::TrackingGeometry>
  buildGenericDetector(Acts::Logging::Level surfaceLLevel = Acts::Logging::INFO,
                       Acts::Logging::Level layerLLevel   = Acts::Logging::INFO,
                       Acts::Logging::Level volumeLLevel  = Acts::Logging::INFO,
                       size_t               version       = 0);

  /// Helper method for positioning
  /// @param radius is the cylinder radius
  /// @param zStagger is the radial staggering along z
  /// @param moduleHalfLength is the module length (longitudinal)
  /// @param lOverlap is the overlap of the modules (longitudinal)
  /// @binningSchema is the way the bins are laid out rphi x z
  std::vector<Acts::Vector3D>
  modulePositionsCylinder(double radius,
                          double zStagger,
                          double moduleHalfLength,
                          double lOverlap,
                          const std::pair<int, int>& binningSchema);

  /// Helper method for positioning
  /// @param z is the z position of the ring
  /// @param radius is the ring radius
  /// @param phiStagger is the radial staggering along phi
  /// @param lOverlap is the overlap of the modules
  /// @parm nPhiBins is the number of bins in phi
  std::vector<Acts::Vector3D>
  modulePositionsRing(double z,
                      double radius,
                      double phiStagger,
                      double phiSubStagger,
                      int    nPhiBins);

  /// Helper method for positioning
  /// @param z is the nominal z posiiton of the dis
  /// @param ringStagger is the staggering of the different rings
  /// @param phiStagger is the staggering on a ring in phi : it is even/odd
  /// @param phiSubStagger is the sub staggering on a ring in phi : it affects
  /// 0/4/8 and 3/6
  /// @param innerRadius is the inner Radius for the disc
  /// @param outerRadius is the outer Radius for the disc
  /// @param discBinning is the binning setup in r, phi
  /// @param moduleHalfLength is pair of phibins and module length
  std::vector<std::vector<Acts::Vector3D>>
  modulePositionsDisc(double                     z,
                      double                     ringStagger,
                      std::vector<double>        phiStagger,
                      std::vector<double>        phiSubStagger,
                      double                     innerRadius,
                      double                     outerRadius,
                      const std::vector<size_t>& discBinning,
                      const std::vector<double>& moduleHalfLength);

}  // end of namespace Generic
}  // end of namespace FW
