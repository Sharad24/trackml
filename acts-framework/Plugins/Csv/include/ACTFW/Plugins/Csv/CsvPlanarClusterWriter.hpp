// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Plugins/Digitization/PlanarModuleCluster.hpp>
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WriterT.hpp"

namespace FW {
namespace Csv {

  /// Write out a planar cluster collection in comma-separated-value format.
  ///
  /// This writes one file per event into the configured output directory. By
  /// default it writes to the current working directory. Files are named
  /// using the following schema
  ///
  ///     event000000001-hits.csv
  ///     event000000002-hits.csv
  ///
  /// and each line in the file corresponds to one hit/cluster.
  class CsvPlanarClusterWriter
      : public WriterT<DetectorData<geo_id_value, Acts::PlanarModuleCluster>>
  {
  public:
    using Base = WriterT<DetectorData<geo_id_value, Acts::PlanarModuleCluster>>;

    struct Config
    {
      std::string collection;           ///< which collection to write
      std::string outputDir;            ///< where to place output files
      size_t      outputPrecision = 6;  ///< floating point precision
    };

    /// Constructor with
    /// @param cfg configuration struct
    /// @param output logging level
    CsvPlanarClusterWriter(const Config&        cfg,
                           Acts::Logging::Level level = Acts::Logging::INFO);

  protected:
    /// This implementation holds the actual writing method
    /// and is called by the WriterT<>::write interface
    ProcessCode
    writeT(const AlgorithmContext& ctx,
           const DetectorData<geo_id_value, Acts::PlanarModuleCluster>&
               clusters) final override;

  private:
    Config m_cfg;
  };
}  // namespace Csv
}  // namespace FW
