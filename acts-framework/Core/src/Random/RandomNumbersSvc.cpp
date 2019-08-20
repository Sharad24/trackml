// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
//  RandomNumbersSvc.cpp
//  ACTFW
//
//  Created by Andreas Salzburger on 17/05/16.
//
//

#include "ACTFW/Random/RandomNumbersSvc.hpp"

FW::RandomNumbersSvc::RandomNumbersSvc(
    const Config&                       cfg,
    std::unique_ptr<const Acts::Logger> logger)
  : m_cfg(cfg), m_logger(std::move(logger))
{
}

std::string
FW::RandomNumbersSvc::name() const
{
  return "RandomNumbersSvc";
}

FW::RandomEngine
FW::RandomNumbersSvc::spawnGenerator(const AlgorithmContext& context) const
{
  return RandomEngine(generateSeed(context));
}

const unsigned int
FW::RandomNumbersSvc::generateSeed(const AlgorithmContext& context) const
{
  // use Cantor pairing function to generate a unique generator id from
  // algorithm and event number to get a consistent seed
  // see https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
  const unsigned int k1 = context.algorithmNumber;
  const unsigned int k2 = context.eventNumber;
  const unsigned int id = (k1 + k2) * (k1 + k2 + 1) / 2 + k2;
  return m_cfg.seed + id;
}
