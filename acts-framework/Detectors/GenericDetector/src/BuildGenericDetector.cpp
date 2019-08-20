// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/GenericDetector/BuildGenericDetector.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include "ACTFW/GenericDetector/GenericLayerBuilder.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Tools/CylinderVolumeBuilder.hpp"
#include "Acts/Tools/CylinderVolumeHelper.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Tools/PassiveLayerBuilder.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Tools/TrackingGeometryBuilder.hpp"
#include "Acts/Tools/TrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/Units.hpp"

namespace FW {

namespace Generic {

  std::unique_ptr<const Acts::TrackingGeometry>
  buildGenericDetector(Acts::Logging::Level surfaceLLevel,
                       Acts::Logging::Level layerLLevel,
                       Acts::Logging::Level volumeLLevel,
                       size_t               stage)
  {
    // configure surface array creator
    auto surfaceArrayCreator
        = std::make_shared<const Acts::SurfaceArrayCreator>(
            Acts::getDefaultLogger("SurfaceArrayCreator", surfaceLLevel));
    // configure the layer creator that uses the surface array creator
    Acts::LayerCreator::Config lcConfig;
    lcConfig.surfaceArrayCreator = surfaceArrayCreator;
    auto layerCreator            = std::make_shared<const Acts::LayerCreator>(
        lcConfig, Acts::getDefaultLogger("LayerCreator", layerLLevel));
    // configure the layer array creator
    auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
        Acts::getDefaultLogger("LayerArrayCreator", layerLLevel));
    // tracking volume array creator
    auto tVolumeArrayCreator
        = std::make_shared<const Acts::TrackingVolumeArrayCreator>(
            Acts::getDefaultLogger("TrackingVolumeArrayCreator", volumeLLevel));
    // configure the cylinder volume helper
    Acts::CylinderVolumeHelper::Config cvhConfig;
    cvhConfig.layerArrayCreator          = layerArrayCreator;
    cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
    auto cylinderVolumeHelper
        = std::make_shared<const Acts::CylinderVolumeHelper>(
            cvhConfig,
            Acts::getDefaultLogger("CylinderVolumeHelper", volumeLLevel));
    //-------------------------------------------------------------------------------------
    // list the volume builders
    std::list<std::shared_ptr<const Acts::ITrackingVolumeBuilder>>
        volumeBuilders;

// a hash include for the Generic Detector : a bit ugly but effective
#include "ACTFW/GenericDetector/GenericDetectorML.ipp"

    //-------------------------------------------------------------------------------------
    // create the tracking geometry
    Acts::TrackingGeometryBuilder::Config tgConfig;
    tgConfig.trackingVolumeBuilders = volumeBuilders;
    tgConfig.trackingVolumeHelper   = cylinderVolumeHelper;
    auto cylinderGeometryBuilder
        = std::make_shared<const Acts::TrackingGeometryBuilder>(
            tgConfig,
            Acts::getDefaultLogger("TrackerGeometryBuilder", volumeLLevel));
    // get the geometry
    auto trackingGeometry = cylinderGeometryBuilder->trackingGeometry();
    /// return the tracking geometry
    return trackingGeometry;
  }

  /// helper method for cylinder
  std::vector<Acts::Vector3D>
  modulePositionsCylinder(double radius,
                          double zStagger,
                          double moduleHalfLength,
                          double lOverlap,
                          const std::pair<int, int>& binningSchema)
  {
    int nPhiBins = binningSchema.first;
    int nZbins   = binningSchema.second;
    // prepare the return value
    std::vector<Acts::Vector3D> mPositions;
    mPositions.reserve(nPhiBins * nZbins);
    // prep work
    double phiStep = 2 * M_PI / (nPhiBins);
    double minPhi  = -M_PI + 0.5 * phiStep;
    double zStart  = -0.5 * (nZbins - 1) * (2 * moduleHalfLength - lOverlap);
    double zStep   = 2 * std::abs(zStart) / (nZbins - 1);
    // loop over the bins
    for (size_t zBin = 0; zBin < size_t(nZbins); ++zBin) {
      // prepare z and r
      double moduleZ = zStart + zBin * zStep;
      double moduleR
          = (zBin % 2) ? radius - 0.5 * zStagger : radius + 0.5 * zStagger;
      for (size_t phiBin = 0; phiBin < size_t(nPhiBins); ++phiBin) {
        // calculate the current phi value
        double modulePhi = minPhi + phiBin * phiStep;
        mPositions.push_back(Acts::Vector3D(
            moduleR * cos(modulePhi), moduleR * sin(modulePhi), moduleZ));
      }
    }
    return mPositions;
  }

  /// helper method for disc
  std::vector<std::vector<Acts::Vector3D>>
  modulePositionsDisc(double                     z,
                      double                     ringStagger,
                      std::vector<double>        phiStagger,
                      std::vector<double>        phiSubStagger,
                      double                     innerRadius,
                      double                     outerRadius,
                      const std::vector<size_t>& discBinning,
                      const std::vector<double>& moduleHalfLength)
  {
    // calculate the radii
    std::vector<double> radii;
    // calculate the radial borders
    std::vector<double> radialBoarders;
    // the radial span of the disc
    double deltaR = outerRadius - innerRadius;
    // quick exits
    if (discBinning.size() == 1) {
      radii.push_back(0.5 * (innerRadius + outerRadius));
      radialBoarders = {innerRadius, outerRadius};
    } else {
      double totalLength = 0;
      // sum up the total length
      for (auto& mhlength : moduleHalfLength) totalLength += 2 * mhlength;
      // now calculate the overlap (equal pay)
      double rOverlap = (totalLength - deltaR) / (moduleHalfLength.size() - 1);
      // and now fill the radii and gaps
      double lastR  = innerRadius;
      double lastHl = 0.;
      double lastOl = 0.;
      // remember the radial boarders
      radialBoarders.push_back(innerRadius);
      // now calculate
      for (auto& mhlength : moduleHalfLength) {
        // calculate the radius
        radii.push_back(lastR + lastHl - lastOl + mhlength);
        lastR  = radii[radii.size() - 1];
        lastOl = rOverlap;
        lastHl = mhlength;
        // and register the radial boarder
        radialBoarders.push_back(lastR + 2 * lastHl - 0.5 * lastOl);
      }
    }
    // now prepare the return method
    std::vector<std::vector<Acts::Vector3D>> mPositions;
    for (size_t ir = 0; ir < radii.size(); ++ir) {
      // generate the z value
      // convention inner ring is closer to origin : makes sense
      double rz = radii.size() == 1 ? z : (ir % 2 ? z + 0.5 * ringStagger
                                                  : z - 0.5 * ringStagger);
      // fill the ring positions
      double psStagger = phiSubStagger.size() ? phiSubStagger[ir] : 0.;
      mPositions.push_back(modulePositionsRing(
          rz, radii[ir], phiStagger[ir], psStagger, discBinning[ir]));
    }
    return mPositions;
  }

  /// Helper method for positioning
  std::vector<Acts::Vector3D>
  modulePositionsRing(double z,
                      double radius,
                      double phiStagger,
                      double phiSubStagger,
                      int    nPhiBins)
  {
    // create and fill the positions
    std::vector<Acts::Vector3D> rPositions;
    rPositions.reserve(nPhiBins);
    // prep work
    double phiStep = 2 * M_PI / (nPhiBins);
    double minPhi  = -M_PI + 0.5 * phiStep;
    // phi loop
    for (size_t iphi = 0; iphi < size_t(nPhiBins); ++iphi) {
      // if we have a phi sub stagger presents
      double rzs = 0.;
      // phi stagger affects 0 vs 1, 2 vs 3 ... etc
      // -> only works if it is a %4
      // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
      if (phiSubStagger != 0. && !(nPhiBins % 4)) {
        // switch sides
        if (!(iphi % 4)) {
          rzs = phiSubStagger;
        } else if (!((iphi + 1) % 4)) {
          rzs = -phiSubStagger;
        }
      }
      // the module phi
      double phi = minPhi + iphi * phiStep;
      // main z position depending on phi bin
      double rz = iphi % 2 ? z - 0.5 * phiStagger : z + 0.5 * phiStagger;
      rPositions.push_back(
          Acts::Vector3D(radius * cos(phi), radius * sin(phi), rz + rzs));
    }
    return rPositions;
  }

}  // end of namespace Generic

}  // end of namespace FW
