// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-10-26 Initial version
/// @author Hadrien Grasland
/// @author Moritz Kiehn <msmk@cern.ch>

#ifndef ACTFW_ALGORITHMCONTEXT_H
#define ACTFW_ALGORITHMCONTEXT_H

#include <memory>

namespace FW {

class WhiteBoard;

/// Aggregated information to run one algorithm over one event.
struct AlgorithmContext
{
  size_t      algorithmNumber;  ///< Unique algorithm identifier
  size_t      eventNumber;      ///< Unique event identifier
  WhiteBoard& eventStore;       ///< Per-event data store
};

}  // namespace FW

#endif  // ACTFW_ALGORITHMCONTEXT_H
