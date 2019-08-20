// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-05-23
/// @author Andreas Salburger
/// @author Moritz Kiehnn <msmk@cern.ch>

#ifndef ACTFW_ISERVICE_H
#define ACTFW_ISERVICE_H

#include <string>

#include "ACTFW/Framework/ProcessCode.hpp"

namespace FW {

/// Interface for common services.
///
/// @warning Use of this interface is **DEPRECATED**. Its only remaining use
///          case is to support legacy writers based on the IWriterT subclass,
///          which should eventually be turned into descendents of the
///          IWriter/WriterT classes. At this point, IService will be removed.
///
/// @todo Remove once all writers have been migrated to IWriter and WriterT
///
class IService
{
public:
  /// Virtual Destructor
  virtual ~IService() = default;

  /// Framework name() method
  virtual std::string
  name() const = 0;

  /// Interface hook to run some code to be executed once all events of a job
  /// have been processed, typically used for commiting writes to a file
  virtual ProcessCode
  endRun()
  {
    return ProcessCode::SUCCESS;
  }
};

}  // namespace FW

#endif  // ACTFW_ISERVICE_H
