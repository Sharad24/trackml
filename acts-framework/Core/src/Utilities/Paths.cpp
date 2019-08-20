// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Utilities/Paths.hpp"
#include <cstdio>

std::string
FW::joinPaths(const std::string& dir, const std::string& name)
{
  if (dir.empty()) {
    return name;
  } else {
    return dir + '/' + name;
  }
}

std::string
FW::perEventFilepath(const std::string& dir,
                     const std::string& name,
                     size_t             event)
{
  char prefix[1024];

  snprintf(prefix, sizeof(prefix), "event%09zu-", event);

  if (dir.empty()) {
    return prefix + name;
  } else {
    return dir + '/' + prefix + name;
  }
}
