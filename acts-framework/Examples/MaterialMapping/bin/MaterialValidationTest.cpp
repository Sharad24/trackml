// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "../src/MaterialValidation.hpp"
#include "ACTFW/Extrapolation/ExtrapolationUtils.hpp"
#include "ACTFW/Framework/Sequencer.hpp"

int
main()
{
  // set geometry building logging level
  Acts::Logging::Level geo - surface - loglevel = Acts::Logging::INFO;
  Acts::Logging::Level layerLogLevel            = Acts::Logging::INFO;
  Acts::Logging::Level volumeLogLevel           = Acts::Logging::INFO;
  // set extrapolation logging level
  Acts::Logging::Level eLogLevel = Acts::Logging::INFO;

  // create the tracking geometry as a shared pointer
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry
      = FW::Generic::buildGenericDetector(
          geo - surface - loglevel, layerLogLevel, volumeLogLevel, 3);

  // set up the magnetic field
  std::shared_ptr<Acts::ConstantBField> magField(
      new Acts::ConstantBField{{0., 0., 2. * Acts::units::_T}});

  // EXTRAPOLATOR - set up the extrapolator
  std::shared_ptr<Acts::IExtrapolationEngine> extrapolationEngine
      = FWE::initExtrapolator(tGeometry, magField, eLogLevel);
}
