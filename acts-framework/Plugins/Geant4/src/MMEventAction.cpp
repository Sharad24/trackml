// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/MMEventAction.hpp"
#include <stdexcept>
#include "ACTFW/Plugins/Geant4/MMPrimaryGeneratorAction.hpp"
#include "ACTFW/Plugins/Geant4/MMSteppingAction.hpp"
#include "G4Event.hh"
#include "G4RunManager.hh"

FW::G4::MMEventAction* FW::G4::MMEventAction::fgInstance = nullptr;

FW::G4::MMEventAction*
FW::G4::MMEventAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

FW::G4::MMEventAction::MMEventAction() : G4UserEventAction()
{
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::G4::MMEventAction::~MMEventAction()
{
  fgInstance = nullptr;
}

void
FW::G4::MMEventAction::BeginOfEventAction(const G4Event*)
{
  // reset the collection of material steps
  MMSteppingAction::Instance()->Reset();
}

void
FW::G4::MMEventAction::EndOfEventAction(const G4Event* event)
{
  const auto*          rawPos = event->GetPrimaryVertex();
  const Acts::Vector3D pos(rawPos->GetX0(), rawPos->GetY0(), rawPos->GetZ0());
  // access the initial direction of the track
  G4ThreeVector rawDir = MMPrimaryGeneratorAction::Instance()->direction();
  const Acts::Vector3D dir(rawDir.x(), rawDir.y(), rawDir.z());
  // create the RecordedMaterialTrack
  Acts::RecordedMaterialTrack mtrecord(
      pos, dir, MMSteppingAction::Instance()->materialSteps());
  // write out the RecordedMaterialTrack of one event
  m_records.push_back(mtrecord);
}

void
FW::G4::MMEventAction::Reset()
{
}
