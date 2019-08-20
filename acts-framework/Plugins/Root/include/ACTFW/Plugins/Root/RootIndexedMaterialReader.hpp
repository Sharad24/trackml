// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_PLUGINS_ROOT_INDEXEDMATERIALREADER_H
#define ACTFW_PLUGINS_ROOT_INDEXEDMATERIALREADER_H

#include <mutex>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Readers/IReaderT.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

class TFile;

namespace FW {

namespace Root {

  /// @class RootMaterialTrackReader
  ///
  /// @brief Reads in MaterialTrack entities from a root file
  ///
  /// This service is the root implementation of the ImaterialTrackReader.
  /// It reads in a vector of MaterialTrack entities from a given root tree
  /// of a given root file. The input file and tree are set over the
  /// configuration
  /// object.
  class RootIndexedMaterialReader
      : public FW::IReaderT<Acts::IndexedSurfaceMaterial>
  {
  public:
    /// @class Config
    /// Configuration of the Reader
    class Config
    {
    public:
      /// The name of the input tree
      std::string folderNameBase = "Material";
      /// The name of the input file
      std::string fileName;
      /// The default logger
      std::shared_ptr<const Acts::Logger> logger;
      /// The name of the service
      std::string name;

      Config(const std::string&   lname = "MaterialReader",
             Acts::Logging::Level lvl   = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname)
      {
      }
    };

    /// Constructor
    RootIndexedMaterialReader(const Config& cfg);

    /// Virtual destructor
    ~RootIndexedMaterialReader() override;

    /// Framework name() method
    std::string
    name() const final override;

    // clang-format off
  /// @copydoc FW::IReaderT::read(std::vector<Acts::ParticleProperties>&,size_t,const FW::AlgorithmContext*)
    // clang-format on
    FW::ProcessCode
    read(Acts::IndexedSurfaceMaterial& mtrc,
         size_t                        skip    = 0,
         const FW::AlgorithmContext*   context = nullptr) final override;

  private:
    /// The config class
    Config m_cfg;

    /// The input file
    TFile* m_inputFile;

    /// mutex used to protect multi-threaded reads
    std::mutex m_read_mutex;

    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }
  };

  inline std::string
  RootIndexedMaterialReader::name() const
  {
    return m_cfg.name;
  }

}  // namespace Root
}  // namespace FW

#endif  // ACTFW_PLUGINS_ROOT_INDEXEDMATERIALREADER_H
