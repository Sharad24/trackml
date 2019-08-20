// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-05-23 Initial version
/// @date 2017-08-07 Rewrite with new interfaces

#ifndef ACTFW_JSONSPACEPOINTWRITER_H
#define ACTFW_JSONSPACEPOINTWRITER_H

#include <fstream>

#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Utilities/Paths.hpp"

namespace FW {
namespace Json {

  /// Write out a space point collection in JSON format.
  ///
  /// This writes one file per event into the configured output directory. By
  /// default it writes to the current working directory. Files are named
  /// using the following schema
  ///
  ///     event000000001-spacepoints.json
  ///     event000000002-spacepoints.json
  ///
  template <class T>
  class JsonSpacePointWriter : public WriterT<DetectorData<geo_id_value, T>>
  {
  public:
    using Base = WriterT<DetectorData<geo_id_value, T>>;

    struct Config
    {
      std::string collection;           ///< which collection to write
      std::string outputDir;            ///< where to place output files
      size_t      outputPrecision = 6;  ///< floating point precision
    };

    JsonSpacePointWriter(const Config&        cfg,
                         Acts::Logging::Level level = Acts::Logging::INFO);

  protected:
    FW::ProcessCode
    writeT(const FW::AlgorithmContext& ctx,
           const DetectorData<geo_id_value, T>& spacePoints) final override;

  private:
    Config m_cfg;
    // required for C++ to find `logger()` with the default look-up
    const Acts::Logger&
    logger() const
    {
      return Base::logger();
    }
  };

}  // namespace Json
}  // namespace FW

template <class T>
FW::Json::JsonSpacePointWriter<T>::JsonSpacePointWriter(
    const FW::Json::JsonSpacePointWriter<T>::Config& cfg,
    Acts::Logging::Level                             level)
  : Base(cfg.collection, "JsonSpacePointWriter", level), m_cfg(cfg)
{
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
}

template <class T>
FW::ProcessCode
FW::Json::JsonSpacePointWriter<T>::writeT(
    const FW::AlgorithmContext& ctx,
    const DetectorData<geo_id_value, T>& spacePoints)
{
  // open per-event file
  std::string path
      = perEventFilepath(m_cfg.outputDir, "spacepoints.json", ctx.eventNumber);
  std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }

  os << std::setprecision(m_cfg.outputPrecision);
  os << "{\n";

  bool firstVolume = true;
  for (auto& volumeData : spacePoints) {
    geo_id_value volumeID = volumeData.first;

    if (!firstVolume) os << ",\n";
    os << "  \"SpacePoints_" << volumeID << "\" : [\n";

    bool firstPoint = true;
    for (auto& layerData : volumeData.second) {
      for (auto& moduleData : layerData.second) {
        for (auto& data : moduleData.second) {
          // set the comma correctly
          if (!firstPoint) os << ",\n";
          // write the space point
          os << "    [" << data.x() << ", " << data.y() << ", " << data.z()
             << "]";
          firstPoint = false;
        }
      }
    }
    os << "]";
    firstVolume = false;
  }
  os << "\n}\n";

  return ProcessCode::SUCCESS;
}

#endif  // ACTFW_JSONSPACEPOINTWRITER_H
