// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/MMSteppingAction.hpp"
#include <stdexcept>
#include "Acts/Utilities/Units.hpp"
#include "G4Material.hh"
#include "G4Step.hh"

FW::G4::MMSteppingAction* FW::G4::MMSteppingAction::fgInstance = nullptr;

FW::G4::MMSteppingAction*
FW::G4::MMSteppingAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

FW::G4::MMSteppingAction::MMSteppingAction() : G4UserSteppingAction(), m_steps()
// m_volMgr(MaterialRunAction::Instance()->getGeant4VolumeManager())
{
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::G4::MMSteppingAction::~MMSteppingAction()
{
  fgInstance = nullptr;
}

void
FW::G4::MMSteppingAction::UserSteppingAction(const G4Step* step)
{
  // get the material
  G4Material* material = step->GetPreStepPoint()->GetMaterial();

  if (material && material->GetName() != "Vacuum"
      && material->GetName() != "Air") {
    // go through the elements of the material & weigh it with its fraction
    const G4ElementVector* elements  = material->GetElementVector();
    const G4double*        fraction  = material->GetFractionVector();
    size_t                 nElements = material->GetNumberOfElements();
    double                 A         = 0.;
    double                 Z         = 0.;
    double                 X0        = material->GetRadlen();
    double                 L0        = material->GetNuclearInterLength();
    double rho        = material->GetDensity() * CLHEP::mm3 / CLHEP::gram;
    double steplength = step->GetStepLength() / CLHEP::mm;
    if (nElements == 1) {
      A = material->GetA() * CLHEP::mole / CLHEP::gram;
      Z = material->GetZ();
    } else {
      for (size_t i = 0; i < nElements; i++) {
        A += elements->at(i)->GetA() * fraction[i];
        Z += elements->at(i)->GetZ() * fraction[i];
      }
      if (A != 0.) A /= nElements;
      if (A != 0.) Z /= nElements;
    }

    /*   G4cout << *material << G4endl;
       G4cout << "----G4StepMaterial----" << G4endl;
       /// @TODO remove output after testing
       G4cout << "Material: " << material->GetName() << G4endl;
       G4cout << "X0: " << X0 << G4endl;
       G4cout << "L0: " << L0 << G4endl;
       G4cout << "A: " << A << G4endl;
       G4cout << "Z: " << Z << G4endl;
       G4cout << "rho: " << rho << G4endl;
       G4cout << "steplength: " << steplength << G4endl;*/

    // create the RecordedMaterialProperties
    const auto&          rawPos = step->GetPreStepPoint()->GetPosition();
    const Acts::Vector3D pos(rawPos.x(), rawPos.y(), rawPos.z());
    const Acts::MaterialProperties matprop(X0, L0, A, Z, rho, steplength);
    const Acts::RecordedMaterialProperties recMatProp(matprop, pos);
    m_steps.push_back(recMatProp);
  }
}

void
FW::G4::MMSteppingAction::Reset()
{
  m_steps.clear();
}
