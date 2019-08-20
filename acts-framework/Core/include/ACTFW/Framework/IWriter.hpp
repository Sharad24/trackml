// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2017-07-25
/// @author Moritz Kiehnn <msmk@cern.ch>

#ifndef ACTFW_IWRITER_H
#define ACTFW_IWRITER_H

#include <string>

#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// Interface for writing data.
class IWriter
{
public:
  /// Virtual destructor
  virtual ~IWriter() = default;

  /// Provide the name of the writer
  virtual std::string
  name() const = 0;

  /// Finish the run (e.g. aggregate statistics, write down output, close files)
  virtual ProcessCode
  endRun()
      = 0;

  /// write data to the output stream
  virtual ProcessCode
  write(const AlgorithmContext& context)
      = 0;
};

}  // namespace FW

#endif  // ACTFW_IWRITER_H
