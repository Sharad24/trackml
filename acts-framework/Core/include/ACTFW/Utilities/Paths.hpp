// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2017-08-03
/// @author Moritz Kiehnn <msmk@cern.ch>

#ifndef ACTFW_PATHS_H
#define ACTFW_PATHS_H

#include <string>

namespace FW {

/// Join dir and name into one paths with correct handling of empty dirs.
std::string
joinPaths(const std::string& dir, const std::string& name);

/// Construct a file path of the form `[<dir>/]event<XXXXX>-<name>`.
///
/// @params dir output directory, unused if empty
/// @params name basic filename
/// @params event event number
std::string
perEventFilepath(const std::string& dir, const std::string& name, size_t event);

}  // namespace FW

#endif  // ACTFW_PATHS_H
