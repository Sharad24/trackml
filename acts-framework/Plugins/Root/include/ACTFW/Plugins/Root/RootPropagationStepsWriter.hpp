// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include "ACTFW/Framework/WriterT.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "TFile.h"
#include "TTree.h"

namespace FW {

namespace Root {

  using PropagationSteps = std::vector<Acts::detail::Step>;

  /// @class RootPropagationStepsWriter
  ///
  /// Write out the steps of test propgations for stepping validation,
  /// each step sequence is one entry in the  in the root file for optimised
  /// data writing speed.
  /// The event number is part of the written data.
  ///
  /// A common file can be provided for to the writer to attach his TTree,
  /// this is done by setting the Config::rootFile pointer to an existing file
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  class RootPropagationStepsWriter
      : public WriterT<std::vector<PropagationSteps>>
  {
  public:
    using Base = WriterT<std::vector<PropagationSteps>>;

    /// @brief The nested configuration struct
    struct Config
    {
      std::string collection
          = "propagation_steps";          ///< particle collection to write
      std::string filePath = "";          ///< path of the output file
      std::string fileMode = "RECREATE";  ///< file access mode
      std::string treeName = "propagation_steps";  ///< name of the output tree
      TFile*      rootFile = nullptr;              ///< common root file
    };

    /// Constructor with
    /// @param cfg configuration struct
    /// @param output logging level
    RootPropagationStepsWriter(const Config&        cfg,
                               Acts::Logging::Level level
                               = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootPropagationStepsWriter() override;

    /// End-of-run hook
    ProcessCode
    endRun() final override;

  protected:
    /// This implementation holds the actual writing method
    /// and is called by the WriterT<>::write interface
    ///
    /// @param ctx The Algorithm context with per event information
    /// @param steps is the data to be written out
    ProcessCode
    writeT(const AlgorithmContext&              ctx,
           const std::vector<PropagationSteps>& steps) final override;

  private:
    Config             m_cfg;          ///< the configuration object
    std::mutex         m_writeMutex;   ///< protect multi-threaded writes
    TFile*             m_outputFile;   ///< the output file name
    TTree*             m_outputTree;   ///< the output tree
    int                m_eventNr;      ///< the event number of
    std::vector<int>   m_volumeID;     ///< volume identifier
    std::vector<int>   m_boundaryID;   ///< boundary identifier
    std::vector<int>   m_layerID;      ///< layer identifier if
    std::vector<int>   m_approachID;   ///< surface identifier
    std::vector<int>   m_sensitiveID;  ///< surface identifier
    std::vector<float> m_x;            ///< global x
    std::vector<float> m_y;            ///< global y
    std::vector<float> m_z;            ///< global z
    std::vector<float> m_dx;           ///< global direction x
    std::vector<float> m_dy;           ///< global direction y
    std::vector<float> m_dz;           ///< global direction z
    std::vector<int>   m_step_type;    ///< step type
    std::vector<float> m_step_acc;     ///< accuracy
    std::vector<float> m_step_act;     ///< actor check
    std::vector<float> m_step_abt;     ///< aborter
    std::vector<float> m_step_usr;     ///< user
  };

}  // namespace Root
}  // namespace FW
