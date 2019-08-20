// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootMaterialTrackWriter.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialStep.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>

#include "TFile.h"

FW::Root::RootMaterialTrackWriter::RootMaterialTrackWriter(
    const FW::Root::RootMaterialTrackWriter::Config& cfg)
  : FW::IWriterT<Acts::MaterialTrack>()
  , m_cfg(cfg)
  , m_outputFile(nullptr)
  , m_outputTree(nullptr)
  , m_trackRecord()
{
  // An input collection name and tree name must be specified
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  } else if (m_cfg.fileName.empty()) {
    throw std::invalid_argument("Missing file name");
  } else if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  } else if (m_cfg.name.empty()) {
    throw std::invalid_argument("Missing service name");
  }

  // Setup ROOT I/O
  m_outputFile = TFile::Open(m_cfg.fileName.c_str(), "recreate");
  if (!m_outputFile) {
    throw std::ios_base::failure("Could not open '" + m_cfg.fileName);
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (!m_outputTree) throw std::bad_alloc();

  // Create a branch with the MaterialTrack entities
  m_outputTree->Branch("MaterialTrack", &m_trackRecord);
}

FW::Root::RootMaterialTrackWriter::~RootMaterialTrackWriter()
{
  m_outputFile->Close();
}

FW::ProcessCode
FW::Root::RootMaterialTrackWriter::endRun()
{
  // write the tree and close the file
  ACTS_INFO("Writing ROOT output File : " << m_cfg.fileName);
  m_outputFile->cd();
  m_outputTree->Write();
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Root::RootMaterialTrackWriter::write(const Acts::MaterialTrack& mtrecord)
{

  // lock the mutex
  std::lock_guard<std::mutex> lock(m_write_mutex);

  // setting the parameters
  m_trackRecord = mtrecord;

  // write to
  m_outputTree->Fill();

  // return success
  return FW::ProcessCode::SUCCESS;
}
