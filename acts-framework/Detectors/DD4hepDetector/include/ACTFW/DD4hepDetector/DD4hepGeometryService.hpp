// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/GeometryInterfaces/IDD4hepService.hpp"
#include "ACTFW/GeometryInterfaces/ITGeoService.hpp"
#include "ACTFW/GeometryInterfaces/ITrackingGeometryService.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "TGeoNode.h"

namespace FW {

namespace DD4hep {

  /// @class DD4hepGeometryService
  ///
  /// @brief service creating geometries from dd4hep input
  ///
  /// The DD4hepGeometryService creates the DD4hep, the TGeo and the ACTS
  /// TrackingGeometry
  /// from DD4hep xml input. The geometries are created only on demand.

  class DD4hepGeometryService : public FW::IDD4hepService,
                                public FW::ITGeoService,
                                public FW::ITrackingGeometryService
  {
  public:
    /// @class Config
    /// nested config file of the DD4hepGeometryService
    class Config
    {
    public:
      /// The default logger
      std::shared_ptr<const Acts::Logger> logger;
      /// XML-file with the detector description
      std::vector<std::string> xmlFileNames;
      /// Logger for the geometry transformation
      Acts::Logging::Level lvl;
      /// The name of the service
      std::string name;
      /// Binningtype in phi
      Acts::BinningType bTypePhi;
      /// Binningtype in r
      Acts::BinningType bTypeR;
      /// Binningtype in z
      Acts::BinningType bTypeZ;
      /// The tolerance added to the geometrical extension in r
      /// of the layers contained to build the volume envelope around
      /// @note this parameter only needs to be set if the volumes containing
      /// the
      /// layers (e.g. barrel, endcap volumes) have no specific shape
      /// (assemblies)
      double envelopeR;
      /// The tolerance added to the geometrical extension in z
      /// of the layers contained to build the volume envelope around
      /// @note this parameter only needs to be set if the volumes containing
      /// the
      /// layers (e.g. barrel, endcap volumes) have no specific shape
      /// (assemblies)
      double envelopeZ;

      double defaultLayerThickness;

      std::function<void(std::vector<dd4hep::DetElement>& detectors)>
          sortDetectors;

      Config(const std::string&   lname = "DD4hepGeometryService",
             Acts::Logging::Level level = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, level))
        , xmlFileNames()
        , lvl(level)
        , name(lname)
        , bTypePhi(Acts::equidistant)
        , bTypeR(Acts::equidistant)
        , bTypeZ(Acts::equidistant)
        , envelopeR(0.)
        , envelopeZ(0.)
      {
      }
    };
    /// Constructor
    DD4hepGeometryService(const Config& cfg);

    /// Virtual destructor
    ~DD4hepGeometryService() override;

    /// Framework name() method
    std::string
    name() const final override;

    /// Interface method to access the DD4hep geometry
    /// @return The world DD4hep DetElement
    dd4hep::DetElement
    dd4hepGeometry() final override;

    /// Interface method to Access the TGeo geometry
    /// @return The world TGeoNode (physical volume)
    TGeoNode*
    tgeoGeometry() final override;

    /// Interface method to access to the interface of the DD4hep geometry
    dd4hep::Detector*
    lcdd() final override;

    /// Interface method to access the ACTS TrackingGeometry
    std::unique_ptr<const Acts::TrackingGeometry>
    trackingGeometry() final override;

  private:
    /// Private method to initiate building of the DD4hep geometry
    FW::ProcessCode
    buildDD4hepGeometry();

    /// Private method to initiate building of the ACTS tracking geometry
    FW::ProcessCode
    buildTrackingGeometry();

    /// The config class
    Config m_cfg;
    /// Pointer to the interface to the DD4hep geometry
    dd4hep::Detector* m_lcdd;
    /// The world DD4hep DetElement
    dd4hep::DetElement m_dd4hepGeometry;
    /// The ACTS TrackingGeometry
    std::unique_ptr<const Acts::TrackingGeometry> m_trackingGeometry;

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }
  };

}  // namespace DD4hep
}  // namespace FW
