// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <TFile.h>
#include <TTree.h>
#include <mutex>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

namespace Root {

  /// @class ExtrapolationCellWriter
  ///
  /// A root based implementation to write out extrapolation steps.
  /// This is the Legacy equivalent of the PropgationSteps writer.
  ///
  /// The event number is part of the written data.
  ///
  /// A common file can be provided for to the writer to attach his TTree,
  /// this is done by setting the Config::rootFile pointer to an existing file
  ///
  /// @tparam parameters_t Type of the track parameters
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  template <typename parameters_t>
  class RootExCellWriter
      : public FW::WriterT<std::vector<Acts::ExtrapolationCell<parameters_t>>>
  {
  public:
    using Base
        = FW::WriterT<std::vector<Acts::ExtrapolationCell<parameters_t>>>;

    using FW::WriterT<std::vector<Acts::ExtrapolationCell<parameters_t>>>::
        logger;

    ///  @struct ExtrapolationStep
    ///  this holds the information to be written out
    struct ExtrapolationStep
    {
      float x, y, z;     ///< position (global)
      float px, py, pz;  ///< momentum
      float type;        ///< type of the step
    };

    /// @brief  The nested config class
    struct Config
    {
    public:
      std::string collection;             ///< particle collection to write
      std::string filePath;               ///< path of the output file
      std::string fileMode = "RECREATE";  ///< file access mode
      std::string treeName
          = "extrapolation_cells";       ///< name of the output tree
      TFile* rootFile       = nullptr;   ///< common root file
      bool   writeSensitive = true;      ///< indiciation to write out sensitive
      bool   writeMaterial  = true;      ///< indiciation to write out material
      bool   writePassive   = true;      ///< indiciation to write out passive
      bool   writeBoundary  = true;      ///< indiciation to write out boundary
      unsigned int reservedSteps = 100;  ///< number of steps to be expected
    };

    /// Constructor
    /// @param cfg is the configuration class
    /// @param level The log level of the writer
    RootExCellWriter(const Config&        cfg,
                     Acts::Logging::Level level = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootExCellWriter() override;

    /// End-of-run hook
    ProcessCode
    endRun() final override;

  protected:
    /// The protected writeT method, called by the WriterT base
    ///
    /// @tparam parameters_t Type of the parameters object
    ///
    /// @param [in] ctx is the algorithm context for event consistency
    /// @param [in] ecells are the celss to be written out
    ProcessCode
    writeT(const FW::AlgorithmContext&                               ctx,
           const std::vector<Acts::ExtrapolationCell<parameters_t>>& ecells)
        final override;

    Config             m_cfg;         ///< the config class
    std::mutex         m_writeMutex;  ///< protect multi-threaded writes
    TFile*             m_outputFile{nullptr};  ///< the output file
    TTree*             m_outputTree{nullptr};  ///< the output tree
    int                m_eventNr;              ///< the event number of
    float              m_eta;                  ///< global eta start
    float              m_phi;                  ///< global phi start
    float              m_materialX0;           ///< material in X0
    float              m_materialL0;           ///< material in L0
    std::vector<float> m_s_positionX;   ///< global position x of the step
    std::vector<float> m_s_positionY;   ///< global position y of the step
    std::vector<float> m_s_positionZ;   ///< global position z of the step
    std::vector<float> m_s_positionR;   ///< global position z of the step
    std::vector<float> m_s_materialX0;  ///< step material X0
    std::vector<float> m_s_materialL0;  ///< step material L0
    std::vector<int>   m_s_material;    ///< type of the step: material
    std::vector<int>   m_s_boundary;    ///< type of the step: boundary
    std::vector<int>   m_s_sensitive;   ///< type of the step: sensitive
    std::vector<int>   m_s_volumeID;    ///< volume identification
    std::vector<int>   m_s_layerID;     ///< layer identification
    std::vector<int>   m_s_surfaceID;   ///< surface identification
    std::vector<float>
        m_s_localposition0;  ///< local position - first coordinate
    std::vector<float>
        m_s_localposition1;  ///< local position - second coordinate
    int m_hits;              ///< number of hits in sensitive material
  };

}  // namespace Root
}  // namespace FW

#include "RootExCellWriter.ipp"
