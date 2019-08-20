// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
//  RandomNumbersDistributions.cpp
//  ACTFW
//
//  Created by Hadrien Grasland on 27/06/17.
//
//

#include "ACTFW/Random/RandomNumberDistributions.hpp"

FW::LandauDist::param_type::param_type(double mean, double scale)
  : mean(mean), scale(scale)
{
}

bool
FW::LandauDist::param_type::operator==(const param_type& other) const
{
  return (mean == other.mean) && (scale == other.scale);
}

FW::LandauDist::LandauDist(double mean, double scale) : m_cfg(mean, scale)
{
}

FW::LandauDist::LandauDist(const param_type& cfg) : m_cfg(cfg)
{
}

FW::LandauDist::result_type
FW::LandauDist::min() const
{
  return -std::numeric_limits<double>::infinity();
}

FW::LandauDist::result_type
FW::LandauDist::max() const
{
  return std::numeric_limits<double>::infinity();
}

bool
FW::LandauDist::operator==(const LandauDist& other) const
{
  return (m_cfg == other.m_cfg);
}
