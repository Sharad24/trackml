// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/MMRunAction.hpp"
#include <stdexcept>
#include "ACTFW/Plugins/Geant4/MMEventAction.hpp"
#include "G4Run.hh"

FW::G4::MMRunAction* FW::G4::MMRunAction::fgInstance = nullptr;

FW::G4::MMRunAction::MMRunAction() : G4UserRunAction()
{
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::G4::MMRunAction::~MMRunAction()
{
  fgInstance = nullptr;
}

FW::G4::MMRunAction*
FW::G4::MMRunAction::Instance()
{
  return fgInstance;
}

void
FW::G4::MMRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  // initialize event cumulative quantities
  MMEventAction::Instance()->Reset();
}

void
FW::G4::MMRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Print
  G4cout << "\n--------------------End of Run------------------------------\n"
         << "\n------------------------------------------------------------\n"
         << G4endl;
}
