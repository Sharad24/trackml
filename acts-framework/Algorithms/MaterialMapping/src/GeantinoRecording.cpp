// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/MaterialMapping/GeantinoRecording.hpp"
#include <iostream>
#include <stdexcept>
#include "ACTFW/Plugins/Geant4/MMDetectorConstruction.hpp"
#include "ACTFW/Plugins/Geant4/MMEventAction.hpp"
#include "ACTFW/Plugins/Geant4/MMPrimaryGeneratorAction.hpp"
#include "ACTFW/Plugins/Geant4/MMRunAction.hpp"
#include "ACTFW/Plugins/Geant4/MMSteppingAction.hpp"
#include "FTFP_BERT.hh"

FW::GeantinoRecording::GeantinoRecording(
    const FW::GeantinoRecording::Config& cnf,
    Acts::Logging::Level                 level)
  : FW::BareAlgorithm("GeantinoRecording", level)
  , m_cfg(cnf)
  , m_runManager(std::make_unique<G4RunManager>())
{
  /// Make sure that a writer was provided in the configuration
  if (!m_cfg.materialTrackWriter) {
    throw std::invalid_argument("Missing material track writer");
  }

  /// Check if the geometry should be accessed over the geant4 service
  if (m_cfg.geant4Service) {
    m_runManager->SetUserInitialization(m_cfg.geant4Service->geant4Geometry());
  } else if (!m_cfg.gdmlFile.empty()) {
    /// Access the geometry from the gdml file
    ACTS_INFO(
        "received Geant4 geometry from GDML file: " << m_cfg.gdmlFile.c_str());
    FW::G4::MMDetectorConstruction* detConstruction
        = new FW::G4::MMDetectorConstruction();
    detConstruction->setGdmlInput(m_cfg.gdmlFile.c_str());
    m_runManager->SetUserInitialization(
        detConstruction);  // constructs detector (calls Construct in
                           // Geant4DetectorConstruction)
  } else {
    throw std::invalid_argument("Missing geometry input for Geant4");
  }

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new FW::G4::MMPrimaryGeneratorAction(
      "geantino", 1000., m_cfg.seed1, m_cfg.seed2));
  FW::G4::MMRunAction* runaction = new FW::G4::MMRunAction();
  m_runManager->SetUserAction(runaction);
  m_runManager->SetUserAction(new FW::G4::MMEventAction());
  m_runManager->SetUserAction(new FW::G4::MMSteppingAction());
  m_runManager->Initialize();
}

FW::ProcessCode FW::GeantinoRecording::execute(FW::AlgorithmContext) const
{

  /// Begin with the simulation
  m_runManager->BeamOn(m_cfg.tracksPerEvent);
  ///
  std::vector<Acts::MaterialTrack> mtrecords
      = FW::G4::MMEventAction::Instance()->MaterialTracks();
  ACTS_INFO("Received " << mtrecords.size()
                        << " MaterialTracks. Writing them now onto file...");
  // write to the file
  for (auto& record : mtrecords) {
    m_cfg.materialTrackWriter->write(record);
  }

  return FW::ProcessCode::SUCCESS;
}
