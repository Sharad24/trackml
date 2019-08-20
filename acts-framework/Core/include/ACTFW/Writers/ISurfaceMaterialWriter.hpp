// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ISurfaceMaterialWriter.h
///////////////////////////////////////////////////////////////////

#ifndef ACTFW_WRITERS_IMATERIALWRITER_H
#define ACTFW_WRITERS_IMATERIALWRITER_H

#include <memory>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace Acts {
class SurfaceMaterial;
class GeometryID;
}

namespace FW {

/// @class ISurfaceMaterialWriter
///
/// Interface class for writing out the material
///

class ISurfaceMaterialWriter : public IService
{
public:
  /// Virtual destructor
  virtual ~ISurfaceMaterialWriter() = default;

  /// Writes out the material map of the layer
  virtual FW::ProcessCode
  write(const Acts::SurfaceMaterial& material,
        const Acts::GeometryID&      geoID,
        const std::string&           name)
      = 0;
};
}
#endif  // ACTFW_WRITERS_IMATERIALWRITER_H
