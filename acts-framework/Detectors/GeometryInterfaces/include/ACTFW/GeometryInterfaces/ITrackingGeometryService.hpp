// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ITrackingGeometryService.hpp
///////////////////////////////////////////////////////////////////

#pragma

#include <memory>
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace Acts {
class TrackingGeometry;
}

namespace FW {

/// @class ITrackingGeometryService
///
/// The ITrackingGeometryService is the interface to access the ACTS tracking
/// geometry.
///
/// @TODO solve problem with double inheritance

class ITrackingGeometryService
{  //: public IService {

public:
  /// virtual destructor
  virtual ~ITrackingGeometryService() = default;

  /// Access to the ACTS TrackingGeometry
  virtual std::unique_ptr<const Acts::TrackingGeometry>
  trackingGeometry() = 0;
};
}
