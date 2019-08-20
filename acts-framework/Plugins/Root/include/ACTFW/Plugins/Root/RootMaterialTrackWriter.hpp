// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_PLUGINS_ROOT_RootMaterialTrackWriter_H
#define ACTFW_PLUGINS_ROOT_RootMaterialTrackWriter_H

#include <mutex>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Writers/IWriterT.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialTrack.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "TTree.h"

namespace FW {

namespace Root {

  /// @class RootMaterialTrackWriter
  ///
  /// @brief Writes out MaterialTrack entities from a root file
  ///
  /// This service is the root implementation of the IWriterT.
  /// It writes out a MaterialTrack which is usually generated from
  /// Geant4 material mapping

  class RootMaterialTrackWriter : public FW::IWriterT<Acts::MaterialTrack>
  {

  public:
    /// @class Config
    /// Configuration of the Writer
    class Config
    {
    public:
      /// The name of the output tree
      std::string treeName;
      /// The name of the output file
      std::string fileName;
      /// The default logger
      std::shared_ptr<const Acts::Logger> logger;
      /// The name of the service
      std::string name;

      Config(const std::string&   lname = "MaterialWriter",
             Acts::Logging::Level lvl   = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname)
      {
      }
    };

    /// Constructor
    RootMaterialTrackWriter(const Config& cfg);

    /// Virtual destructor
    ~RootMaterialTrackWriter() override;

    /// Framework name() method
    std::string
    name() const final override;

    /// Framework intialize method
    FW::ProcessCode
    endRun() final override;

    /// Interface method which writes out the MaterialTrack entities
    /// @param mtrecord is the material track record to be written out
    FW::ProcessCode
    write(const Acts::MaterialTrack& mtrecord) final override;

  private:
    /// The config class
    Config m_cfg;
    /// mutex used to protect multi-threaded writes
    std::mutex m_write_mutex;
    /// The output file name
    TFile* m_outputFile;
    /// The output tree name
    TTree* m_outputTree;
    // The MaterialTrack to be written out
    Acts::MaterialTrack m_trackRecord;

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }
  };

  inline std::string
  RootMaterialTrackWriter::name() const
  {
    return m_cfg.name;
  }

}  // namespace Root
}  // namespace FW

#endif  // ACTFW_PLUGINS_ROOT_RootMaterialTrackWriter_H
