// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootIndexedMaterialWriter.hpp"
#include <ios>
#include <iostream>
#include <stdexcept>
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "TFile.h"
#include "TH2F.h"

FW::Root::RootIndexedMaterialWriter::RootIndexedMaterialWriter(
    const FW::Root::RootIndexedMaterialWriter::Config& cfg)
  : FW::IWriterT<Acts::IndexedSurfaceMaterial>()
  , m_cfg(cfg)
  , m_outputFile(nullptr)
{
  // Validate the configuration
  if (m_cfg.folderNameBase.empty()) {
    throw std::invalid_argument("Missing folder name base");
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
}

FW::Root::RootIndexedMaterialWriter::~RootIndexedMaterialWriter()
{
  m_outputFile->Close();
}

FW::ProcessCode
FW::Root::RootIndexedMaterialWriter::write(
    const Acts::IndexedSurfaceMaterial& ism)
{
  // lock the mutex
  std::lock_guard<std::mutex> lock(m_write_mutex);

  // get the geometry ID
  Acts::GeometryID geoID = ism.first;
  // decode the geometryID
  geo_id_value gvolID = geoID.value(Acts::GeometryID::volume_mask);
  geo_id_value glayID = geoID.value(Acts::GeometryID::layer_mask);
  geo_id_value gappID = geoID.value(Acts::GeometryID::approach_mask);
  geo_id_value gsenID = geoID.value(Acts::GeometryID::sensitive_mask);
  // create the directory
  std::string tdName = m_cfg.folderNameBase.c_str();
  tdName += "_vol" + std::to_string(gvolID);
  tdName += "_lay" + std::to_string(glayID);
  tdName += "_app" + std::to_string(gappID);
  tdName += "_sen" + std::to_string(gsenID);
  // create a new directory
  m_outputFile->mkdir(tdName.c_str());
  m_outputFile->cd(tdName.c_str());

  ACTS_INFO("Writing out map at " << tdName);

  size_t bins0 = 0, bins1 = 0;
  // understand what sort of material you have in mind
  Acts::BinnedSurfaceMaterial* bsm
      = dynamic_cast<Acts::BinnedSurfaceMaterial*>(ism.second);
  if (bsm) {
    // overwrite the bin numbers
    bins0 = bsm->binUtility().bins(0);
    bins1 = bsm->binUtility().bins(1);
  }
  TH2F* t = new TH2F("t",
                     "thickness [mm] ;b0 ;b1",
                     bins0,
                     -0.5,
                     bins0 - 0.5,
                     bins1,
                     -0.5,
                     bins1 - 0.5);
  TH2F* x0 = new TH2F("X0",
                      "X_{0} [mm] ;b0 ;b1",
                      bins0,
                      -0.5,
                      bins0 - 0.5,
                      bins1,
                      -0.5,
                      bins1 - 0.5);
  TH2F* l0 = new TH2F("L0",
                      "#Lambda_{0} [mm] ;b0 ;b1",
                      bins0,
                      -0.5,
                      bins0 - 0.5,
                      bins1,
                      -0.5,
                      bins1 - 0.5);
  TH2F* A = new TH2F("A",
                     "X_{0} [mm] ;b0 ;b1",
                     bins0,
                     -0.5,
                     bins0 - 0.5,
                     bins1,
                     -0.5,
                     bins1 - 0.5);
  TH2F* Z = new TH2F("Z",
                     "#Lambda_{0} [mm] ;b0 ;b1",
                     bins0,
                     -0.5,
                     bins0 - 0.5,
                     bins1,
                     -0.5,
                     bins1 - 0.5);
  TH2F* rho = new TH2F("rho",
                       "#rho [g/mm^3] ;b0 ;b1",
                       bins0,
                       -0.5,
                       bins0 - 0.5,
                       bins1,
                       -0.5,
                       bins1 - 0.5);

  // loop over the material and fill
  for (size_t b0 = 0; b0 < bins0; ++b0)
    for (size_t b1 = 0; b1 < bins1; ++b1) {
      // get the material for the bin
      auto mat = bsm->material(b0, b1);
      if (mat) {
        t->SetBinContent(b0 + 1, b1 + 1, mat->thickness());
        x0->SetBinContent(b0 + 1, b1 + 1, mat->material().X0());
        l0->SetBinContent(b0 + 1, b1 + 1, mat->material().L0());
        A->SetBinContent(b0 + 1, b1 + 1, mat->material().A());
        Z->SetBinContent(b0 + 1, b1 + 1, mat->material().Z());
        rho->SetBinContent(b0 + 1, b1 + 1, mat->material().rho());
      }
    }
  t->Write();
  x0->Write();
  l0->Write();
  A->Write();
  Z->Write();
  rho->Write();

  // return success
  return FW::ProcessCode::SUCCESS;
}
