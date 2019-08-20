// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IMaterialStepWriter.h
///////////////////////////////////////////////////////////////////

#ifndef WRITERS_IMATERIALSTEPWRITER_H
#define WRITERS_IMATERIALSTEPWRITER_H

#include <string>
#include <vector>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
class MaterialStep;
class Surface;
}
namespace FW {

/// @class IMaterialStepWriter
/// Interface to write out a vector of MaterialStep entities per surface
///

class IMaterialStepWriter : public IService
{
public:
  /// Virtual destructor
  virtual ~IMaterialStepWriter() = default;

  /// Writes out the MaterialStep entities
  /// @param surface the underlying Acts::Surface corresponding to the material
  /// steps and layer positions
  /// @param steps all MaterialSteps of this layer/surface
  /// @param realAnsAssignedPos the real and corresponding assigned positions of
  /// the material
  virtual ProcessCode
  write(std::string                           name,
        const Acts::Surface*                  surface,
        const std::vector<Acts::MaterialStep> steps,
        const std::vector<std::pair<const Acts::Vector3D, const Acts::Vector3D>>
            realAndAssignedPos)
      = 0;
};
}
#endif  // WRITERS_IMATERIALSTEPWRITER_H
