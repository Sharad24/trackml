// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/ParametricProcessGenerator.hpp"

#include <random>

FW::ParametricProcessGenerator::ParametricProcessGenerator(
    const FW::ParametricProcessGenerator::Config& cfg)
  : m_cfg(cfg)
{
}

std::vector<FW::Data::SimVertex<FW::Data::SimParticle>>
FW::ParametricProcessGenerator::operator()(FW::RandomEngine& rng) const
{
  using Uniform = std::uniform_real_distribution<double>;

  Uniform d0Dist(m_cfg.d0Range[0], m_cfg.d0Range[1]);
  Uniform z0Dist(m_cfg.z0Range[0], m_cfg.z0Range[1]);
  Uniform phiDist(m_cfg.phiRange[0], m_cfg.phiRange[1]);
  Uniform etaDist(m_cfg.etaRange[0], m_cfg.etaRange[1]);
  Uniform ptDist(m_cfg.ptRange[0], m_cfg.ptRange[1]);
  auto generateChargeSign = [=](RandomEngine& rng) {
    return (m_cfg.randomCharge and (Uniform(0.0, 1.0)(rng) < 0.5)) ? -1.0 : 1.0;
  };

  // create empty process vertex
  Data::SimVertex<Data::SimParticle> process({0.0, 0.0, 0.0});

  for (size_t ip = 0; ip < m_cfg.numParticles; ip++) {
    auto d0    = d0Dist(rng);
    auto z0    = z0Dist(rng);
    auto phi   = phiDist(rng);
    auto eta   = etaDist(rng);
    auto pt    = ptDist(rng);
    auto qsign = generateChargeSign(rng);

    Acts::Vector3D position(d0 * std::sin(phi), d0 * -std::cos(phi), z0);
    Acts::Vector3D momentum(
        pt * std::cos(phi), pt * std::sin(phi), pt * std::sinh(eta));
    process.out.emplace_back(position,
                             momentum,
                             m_cfg.mass,
                             qsign * m_cfg.charge,
                             qsign * m_cfg.pdg);
  }
  return {process};
}
