// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DD4HEP_DD4HEPTOG4SVC_H
#define DD4HEP_DD4HEPTOG4SVC_H

#include <memory>
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/GeometryInterfaces/IDD4hepService.hpp"
#include "ACTFW/GeometryInterfaces/IGeant4Service.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"

namespace FW {

namespace DD4hepG4 {

  /// @class DD4hepToG4Svc
  ///
  /// The DD4hepToG4Svc provides the Geant4 geometry from DD4hep input.
  /// @TODO This class should only be temporary - later it should be replaced by
  /// TGeoToGeant4

  class DD4hepToG4Svc : public FW::IGeant4Service
  {

  public:
    /// @class Config
    /// nested config file of DD4hepToG4Svc
    class Config
    {
    public:
      /// The default logger
      std::shared_ptr<const Acts::Logger> logger;
      /// The Service for accessing the DD4hep geometry
      std::shared_ptr<FW::IDD4hepService> dd4hepService;
      /// The name of the service
      std::string name;

      Config(const std::string&   lname = "MaterialWriter",
             Acts::Logging::Level lvl   = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl))
        , dd4hepService(nullptr)
        , name(lname)
      {
      }
    };
    /// Constructor
    DD4hepToG4Svc(const Config& cfg);

    /// Virtual destructor
    ~DD4hepToG4Svc() override;

    /// Framework name() method
    std::string
    name() const final override;

    /// Interface method to access the geant4 geometry
    /// @return G4VUserDetectorConstruction from which the Geant4 geometry is
    /// constructed
    G4VUserDetectorConstruction*
    geant4Geometry() final override;

  private:
    /// Private method building the detector in Geant4 geometry
    FW::ProcessCode
    buildGeant4Geometry();
    /// The config class
    Config m_cfg;
    /// Pointer to the interface to the DD4hep geometry
    G4VUserDetectorConstruction* m_geant4Geometry;

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }
  };

}  // namespace DD4hepG4
}  // namespace FW

#endif  // DD4HEP_DD4HEPTOG4SVC_H
