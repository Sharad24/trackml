// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <fstream>
#include <iostream>
#include <mutex>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Writers/IWriterT.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {
namespace Csv {

  /// @class CsvSurfaceWriter
  ///
  /// A Cvs surface writer for the geometry
  ///
  class CsvSurfaceWriter : public FW::IWriterT<Acts::Surface>
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
      /// the name of the algorithm
      std::string name;
      /// write sensitive surfaces
      bool outputSensitive = true;
      /// write the layer surface out
      bool outputLayerSurface = false;
      /// write the bounds
      bool outputBounds = true;
      /// output scalor
      double outputScalor = 1.;
      /// precision for out
      unsigned int outputPrecision = 6;
      /// file prefix to be written out
      std::string filePrefix = "";
      /// the output stream
      std::shared_ptr<std::ofstream> outputStream = nullptr;

      /// Constructor of nested config class
      /// @param lname name of the Writer
      /// @param lvl output log level
      Config(const std::string&   lname = "CsvSurfaceWriter",
             Acts::Logging::Level lvl   = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname)
      {
      }
    };

    /// Constructor
    ///
    /// @param cfg is the configuration class
    CsvSurfaceWriter(const Config& cfg);

    /// Framework name() method
    std::string
    name() const final override;

    /// The write interface
    /// @param surface to be written out
    FW::ProcessCode
    write(const Acts::Surface& surface) final override;

    /// write a bit of string
    /// @param is the string to be written
    FW::ProcessCode
    write(const std::string& sinfo);

  private:
    Config     m_cfg;          ///< the config class
    std::mutex m_write_mutex;  ///< mutex to protect multi-threaded writes

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }
  };

  inline FW::ProcessCode
  CsvSurfaceWriter::write(const std::string& sinfo)
  {
    // lock the mutex for writing
    std::lock_guard<std::mutex> lock(m_write_mutex);
    // and write
    (*m_cfg.outputStream) << sinfo;
    return FW::ProcessCode::SUCCESS;
  }

}  // namespace Csv
}  // namespace FW
