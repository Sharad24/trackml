// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ISurfaceMaterialReader.h
///////////////////////////////////////////////////////////////////

#ifndef ACTFW_READERS_IMATERIALREADER_H
#define ACTFW_READERS_IMATERIALREADER_H

#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace Acts {
class SurfaceMaterial;
}

namespace FW {

/// @class ISurfaceMaterialReader
///
/// Interface class for reading in the  material
///

class ISurfaceMaterialReader : public IService
{

public:
  /// Virtual destructor
  virtual ~ISurfaceMaterialReader() = default;

  /// Writes out the material map of the
  virtual ProcessCode
  read()
      = 0;
};
}
#endif  // ACTFW_READERS_IMATERIALREADER_H