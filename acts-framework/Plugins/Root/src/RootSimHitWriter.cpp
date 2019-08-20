// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootSimHitWriter.hpp"
#include <ios>
#include <stdexcept>
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Plugins/Digitization/Segmentation.hpp"

FW::Root::RootSimHitWriter::RootSimHitWriter(
    const FW::Root::RootSimHitWriter::Config& cfg,
    Acts::Logging::Level                      level)
  : Base(cfg.collection, "RootSimHitWriter", level)
  , m_cfg(cfg)
  , m_outputFile(cfg.rootFile)
{
  // An input collection name and tree name must be specified
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_cfg.rootFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (!m_outputFile) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }

  m_outputFile->cd();
  m_outputTree
      = new TTree(m_cfg.treeName.c_str(), "TTree from RootSimHitWriter");
  if (m_outputTree == nullptr) throw std::bad_alloc();

  // Set the branches
  m_outputTree->Branch("volume_id", &m_volumeID);
  m_outputTree->Branch("layer_id", &m_layerID);
  m_outputTree->Branch("surface_id", &m_surfaceID);
  m_outputTree->Branch("g_x", &m_x);
  m_outputTree->Branch("g_y", &m_y);
  m_outputTree->Branch("g_z", &m_z);
  m_outputTree->Branch("d_x", &m_dx);
  m_outputTree->Branch("d_y", &m_dy);
  m_outputTree->Branch("d_z", &m_dz);
  m_outputTree->Branch("value", &m_value);
}

FW::Root::RootSimHitWriter::~RootSimHitWriter()
{
  /// Close the file if it's yours
  if (m_cfg.rootFile == nullptr) {
    m_outputFile->Close();
  }
}

FW::ProcessCode
FW::Root::RootSimHitWriter::endRun()
{
  // Write the tree
  m_outputFile->cd();
  m_outputTree->Write();
  ACTS_VERBOSE("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                           << m_cfg.filePath
                                           << "'");
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Root::RootSimHitWriter::writeT(
    const AlgorithmContext& ctx,
    const FW::DetectorData<geo_id_value, Data::SimHit<Data::SimParticle>>&
        fhits)
{
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Loop over the planar fatras hits in this event
  for (auto& volumeData : fhits) {
    for (auto& layerData : volumeData.second) {
      for (auto& moduleData : layerData.second) {
        for (auto& hit : moduleData.second) {
          // Volume / layer and surface identification
          m_volumeID  = volumeData.first;
          m_layerID   = layerData.first;
          m_surfaceID = moduleData.first;
          m_x         = hit.position.x();
          m_y         = hit.position.y();
          m_z         = hit.position.z();
          m_dx        = hit.direction.x();
          m_dy        = hit.direction.y();
          m_dz        = hit.direction.z();
          m_value     = hit.value;
          // Fill the tree
          m_outputTree->Fill();
        }
      }
    }
  }
  return FW::ProcessCode::SUCCESS;
}
