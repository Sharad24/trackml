// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// Interface for reading data.
class IReader
{
public:
  /// Virtual destructor
  virtual ~IReader() = default;

  /// Provide the name of the reader
  virtual std::string
  name() const = 0;

  /// Skip the next n events.
  virtual ProcessCode
  skip(size_t skip)
      = 0;

  /// Read the next event
  virtual ProcessCode
  read(AlgorithmContext context)
      = 0;

  /// The number of events to be read in
  /// @return number of events
  virtual size_t
  numEvents() const = 0;
};

}  // namespace FW
