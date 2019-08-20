// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MMMaterialStepAction.h
///////////////////////////////////////////////////////////////////

#ifndef ACTFW_PLUGINS_GEANT4_MMSTEPPINGACTION_H
#define ACTFW_PLUGINS_GEANT4_MMSTEPPINGACTION_H

#include <vector>
#include "Acts/Plugins/MaterialMapping/RecordedMaterialTrack.hpp"
#include "G4UserSteppingAction.hh"
#include "globals.hh"

namespace FW {
namespace G4 {

  /// @class MMSteppingAction
  ///
  /// @brief Collects the RecordedMaterialProperties entities
  ///
  /// The MMSteppingAction class is the implementation of the
  /// Geant4 class SteppingAction. It creates extracts the weighted material
  /// of every step and collects all material steps.

  class MMSteppingAction : public G4UserSteppingAction
  {
  public:
    /// Constructor
    MMSteppingAction();

    /// Destructor
    ~MMSteppingAction() override;

    /// Static access method
    static MMSteppingAction*
    Instance();

    /// Interface Method doing the step
    /// @note it creates and collects the RecordedMaterialProperties entities
    /// @param step is the Geant4 step of the particle
    void
    UserSteppingAction(const G4Step* step) final override;

    /// Interface reset method
    /// @note it clears the collected step vector
    void
    Reset();

    /// Access to the collected RecordedMaterialProperties entities
    std::vector<Acts::RecordedMaterialProperties>
    materialSteps()
    {
      return m_steps;
    }

  private:
    /// Instance of the SteppingAction
    static MMSteppingAction* fgInstance;

    /// The collected RecordedMaterialProperties entities
    std::vector<Acts::RecordedMaterialProperties> m_steps;
  };

}  // namespace G4
}  // namespace FW

#endif  // ACTFW_PLUGINS_GEANT4_MMSTEPPINGACTION_H
