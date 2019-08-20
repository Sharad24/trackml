// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "TFile.h"
#include "TTree.h"

namespace FW {

namespace Root {

  /// @class RootSimHitWriter
  ///
  /// Write out a planar cluster collection into a root file
  /// to avoid immense long vectors, each cluster is one entry
  /// in the root file for writing speed optimisation.
  /// The event number is part of the written data.
  ///
  /// A common file can be provided for to the writer to attach his TTree,
  /// this is done by setting the Config::rootFile pointer to an existing file
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  class RootSimHitWriter
      : public WriterT<DetectorData<geo_id_value,
                                    Data::SimHit<Data::SimParticle>>>
  {
  public:
    using Base
        = WriterT<DetectorData<geo_id_value, Data::SimHit<Data::SimParticle>>>;

    /// @brief The nested configuration struct
    struct Config
    {
      std::string collection;             ///< cluster collection to write
      std::string filePath;               ///< path of the output file
      std::string fileMode = "RECREATE";  ///< file access mode
      std::string treeName = "hits";      ///< name of the output tree
      TFile*      rootFile = nullptr;     ///< common root file
    };

    /// Constructor with
    /// @param cfg configuration struct
    /// @param output logging level
    RootSimHitWriter(const Config&        cfg,
                     Acts::Logging::Level level = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootSimHitWriter() override;

    /// End-of-run hook
    ProcessCode
    endRun() final override;

  protected:
    /// This implementation holds the actual writing method
    /// and is called by the WriterT<>::write interface
    ///
    /// @param ctx The Algorithm context with per event information
    /// @param simhits The simulation hits collection to we written out
    ProcessCode
    writeT(const AlgorithmContext& ctx,
           const DetectorData<geo_id_value, Data::SimHit<Data::SimParticle>>&
               simhits) final override;

  private:
    Config     m_cfg;         ///< the configuration object
    std::mutex m_writeMutex;  ///< protect multi-threaded writes
    TFile*     m_outputFile;  ///< the output file
    TTree*     m_outputTree;  ///< the output tree to be written to
    int        m_eventNr;     ///< the event number of
    int        m_volumeID;    ///< volume identifier
    int        m_layerID;     ///< layer identifier
    int        m_surfaceID;   ///< surface identifier
    float      m_x;           ///< global x
    float      m_y;           ///< global y
    float      m_z;           ///< global z
    float      m_dx;          ///< global direction x
    float      m_dy;          ///< global direction y
    float      m_dz;          ///< global direction z
    float      m_value;       ///< value of the hit
  };

}  // namespace Root
}  // namespace FW
