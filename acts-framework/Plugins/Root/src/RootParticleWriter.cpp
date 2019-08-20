// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootParticleWriter.hpp"
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

FW::Root::RootParticleWriter::RootParticleWriter(
    const FW::Root::RootParticleWriter::Config& cfg,
    Acts::Logging::Level                        level)
  : ParticleWriter(cfg.collection, "RootParticleWriter", level)
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
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("eta", &m_eta);
    m_outputTree->Branch("phi", &m_phi);
    m_outputTree->Branch("vx", &m_vx);
    m_outputTree->Branch("vy", &m_vy);
    m_outputTree->Branch("vz", &m_vz);
    m_outputTree->Branch("px", &m_px);
    m_outputTree->Branch("py", &m_py);
    m_outputTree->Branch("pz", &m_pz);
    m_outputTree->Branch("pt", &m_pT);
    m_outputTree->Branch("charge", &m_charge);
    m_outputTree->Branch("mass", &m_mass);
    m_outputTree->Branch("pdg", &m_pdgCode);
    m_outputTree->Branch("barcode", &m_barcode, "barcode/l");
    m_outputTree->Branch("vertex", &m_vertex);
    m_outputTree->Branch("primary", &m_primary);
    m_outputTree->Branch("generation", &m_generation);
    m_outputTree->Branch("secondary", &m_secondary);
    m_outputTree->Branch("process", &m_process);
  }
}

FW::Root::RootParticleWriter::~RootParticleWriter()
{
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode
FW::Root::RootParticleWriter::endRun()
{
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Wrote particles to tree '" << m_cfg.treeName << "' in '"
                                          << m_cfg.filePath
                                          << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Root::RootParticleWriter::writeT(
    const AlgorithmContext&               ctx,
    const std::vector<Data::SimVertex<>>& vertices)
{

  if (m_outputFile == nullptr) return ProcessCode::SUCCESS;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // loop over the process vertices
  for (auto& vertex : vertices) {
    for (auto& particle : vertex.outgoing()) {
      /// collect the information
      m_vx      = particle.position().x();
      m_vy      = particle.position().y();
      m_vz      = particle.position().z();
      m_eta     = eta(particle.momentum());
      m_phi     = phi(particle.momentum());
      m_px      = particle.momentum().x();
      m_py      = particle.momentum().y();
      m_pz      = particle.momentum().z();
      m_pT      = perp(particle.momentum());
      m_charge  = particle.q();
      m_mass    = particle.m();
      m_pdgCode = particle.pdg();

      auto barcode = particle.barcode();
      m_barcode    = barcode;
      // decode using the barcode service
      if (m_cfg.barcodeSvc) {
        // the barcode service
        m_vertex     = m_cfg.barcodeSvc->vertex(barcode);
        m_primary    = m_cfg.barcodeSvc->primary(barcode);
        m_generation = m_cfg.barcodeSvc->generate(barcode);
        m_secondary  = m_cfg.barcodeSvc->secondary(barcode);
        m_process    = m_cfg.barcodeSvc->process(barcode);
      }
      m_outputTree->Fill();
    }
  }

  return ProcessCode::SUCCESS;
}
