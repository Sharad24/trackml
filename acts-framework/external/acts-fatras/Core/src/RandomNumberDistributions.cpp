// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras/Kernel/detail/RandomNumberDistributions.hpp"

Fatras::LandauDist::param_type::param_type(double mean, double scale)
    : mean(mean), scale(scale) {}

bool Fatras::LandauDist::param_type::operator==(const param_type &other) const {
  return (mean == other.mean) && (scale == other.scale);
}

Fatras::LandauDist::LandauDist(double mean, double scale)
    : m_cfg(mean, scale) {}

Fatras::LandauDist::LandauDist(const param_type &cfg) : m_cfg(cfg) {}

Fatras::LandauDist::result_type Fatras::LandauDist::min() const {
  return -std::numeric_limits<double>::infinity();
}

Fatras::LandauDist::result_type Fatras::LandauDist::max() const {
  return std::numeric_limits<double>::infinity();
}

bool Fatras::LandauDist::operator==(const LandauDist &other) const {
  return (m_cfg == other.m_cfg);
}
