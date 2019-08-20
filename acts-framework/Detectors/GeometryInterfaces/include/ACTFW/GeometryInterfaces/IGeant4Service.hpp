// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
/// IGeant4Service.hpp
///////////////////////////////////////////////////////////////////

#pragma once

#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

class G4VUserDetectorConstruction;

namespace FW {

/// @class IGeant4Service
///
/// The IGeant4Service is an interface class to return the Geant4 geometry.

class IGeant4Service : public IService
{

public:
  /// virtual destructor
  virtual ~IGeant4Service() = default;
  /// Access to the geant4 geometry
  /// @return G4VUserDetectorConstruction from which the Geant4 geometry is
  /// constructed
  virtual G4VUserDetectorConstruction*
  geant4Geometry()
      = 0;
};
}
