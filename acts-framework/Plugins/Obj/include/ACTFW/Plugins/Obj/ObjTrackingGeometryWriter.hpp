// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>

#include <fstream>
#include <iostream>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Plugins/Obj/ObjSurfaceWriter.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {
class TrackingVolume;
class TrackingGeometry;
}

namespace FW {

namespace Obj {

  /// @class ObjTrackingGeometryWriter
  ///
  /// An Obj writer for the geometry: TrackingGeometry master
  /// It delegates the writing of surfaces to the surface writers
  class ObjTrackingGeometryWriter : public FW::IWriterT<Acts::TrackingGeometry>
  {
  public:
    // @class Config
    //
    // The nested config class
    class Config
    {
    public:
      /// the default logger
      std::shared_ptr<const Acts::Logger> logger;
      /// the name of the writer
      std::string name = "";
      /// surfaceWriters
      std::vector<std::shared_ptr<ObjSurfaceWriter>> surfaceWriters;
      std::string                                    filePrefix           = "";
      std::string                                    sensitiveGroupPrefix = "";
      std::string                                    layerPrefix          = "";

      Config(const std::string&   lname = "ObjTrackingGeometryWriter",
             Acts::Logging::Level lvl   = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl))
        , name(lname)
        , surfaceWriters()
      {
      }
    };

    /// Constructor
    /// @param cfg is the configuration class
    ObjTrackingGeometryWriter(const Config& cfg);

    /// Framework name() method
    /// @return the name of the tool
    std::string
    name() const final override;

    /// The write interface
    /// @param tGeometry is the geometry to be written out
    /// @return ProcessCode to indicate success/failure
    FW::ProcessCode
    write(const Acts::TrackingGeometry& tGeometry) final override;

  private:
    Config m_cfg;  ///< the config class

    /// process this volume
    /// @param tVolume the volume to be processed
    void
    write(const Acts::TrackingVolume& tVolume);

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }
  };

}  // namespace Obj
}  // namepace FW
