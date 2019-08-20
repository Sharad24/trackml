// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MMEventAction.h
///////////////////////////////////////////////////////////////////

#ifndef ACTFW_PLUGINS_GEANT4_MMEVENTACTION_H
#define ACTFW_PLUGINS_GEANT4_MMEVENTACTION_H

#include <memory>
#include "Acts/Plugins/MaterialMapping/RecordedMaterialTrack.hpp"
#include "G4UserEventAction.hh"
#include "globals.hh"

/// @namespace FW::G4:: Namespace for geant4 material mapping
namespace FW {
namespace G4 {

  class MMSteppingAction;

  /// @class MMEventAction
  ///
  /// @brief Writes out material track records
  ///
  /// The MMEventAction class is the realization of the Geant4 class
  /// G4UserEventAction and is writing out the collected RecordedMaterialTrack
  /// entities needed for material mapping once per event.
  ///
  class MMEventAction : public G4UserEventAction
  {
  public:
    /// Constructor
    MMEventAction();

    /// Virtual destructor
    ~MMEventAction() override;

    /// Static access method
    static MMEventAction*
    Instance();

    /// Interface method for begin of the event
    /// @param event is the G4Event to be processed
    /// @note resets the material step action
    void
    BeginOfEventAction(const G4Event* event) final override;

    /// Interface method for end of event
    /// @param event is the G4Event to be processed
    /// @note this method is writing out the material track records
    void
    EndOfEventAction(const G4Event* event) final override;

    /// Interface method
    /// @note does nothing
    void
    Reset();

    // Access the material track records
    std::vector<Acts::RecordedMaterialTrack> const
    MaterialTracks();

  private:
    /// Instance of the EventAction
    static MMEventAction* fgInstance;

    /// The materialTrackWriter
    std::vector<Acts::RecordedMaterialTrack> m_records;
  };

  inline std::vector<Acts::RecordedMaterialTrack> const
  MMEventAction::MaterialTracks()
  {
    auto rrecords = m_records;
    m_records.clear();
    return rrecords;
  }
}  // namespace G4
}  // namespace FW

#endif  // ACTFW_PLUGINS_GEANT4_MMEVENTACTION_H
