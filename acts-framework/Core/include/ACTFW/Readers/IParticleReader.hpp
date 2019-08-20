// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//  IExtrapolationCellWriter.h
//  ACTS-Development
//
//  Created by Andreas Salzburger on 23/05/16.
//
//

#ifndef ACTFW_READERS_IPARTICLEREADER_H
#define ACTFW_READERS_IPARTICLEREADER_H 1

#include "Acts/EventData/ParticleDefinitions.hpp"

namespace FW {

/// @class IParticleReader
///
/// Interface class that fills a vector of particle
/// proerties for feeding into the fast simulation
///
class IParticleReader
{
public:
  /// Virtual Destructor
  virtual ~IParticleReader() = default;

  /// returns the list of particles to be processed
  /// for the fast simulation
  virtual std::vector<Acts::ParticleProperties>
  particles() const = 0;

  /// skip some events
  virtual void
  skip(size_t nEvents = 1) const = 0;
};
}

#endif  // ACTFW_READERS_IPARTICLEREADER_H
