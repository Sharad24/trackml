// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/ParticleSelector.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FW::ParticleSelector::ParticleSelector(const Config&        cfg,
                                       Acts::Logging::Level level)
  : FW::BareAlgorithm("Selector", level), m_cfg(cfg)
{
  if (m_cfg.input.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
  if (m_cfg.output.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
}

FW::ProcessCode
FW::ParticleSelector::execute(FW::AlgorithmContext ctx) const
{
  const std::vector<Data::SimVertex<Data::SimParticle>>* all = nullptr;
  std::vector<Data::SimVertex<Data::SimParticle>>        selected;

  // get input particles
  if (ctx.eventStore.get(m_cfg.input, all) == ProcessCode::ABORT)
    return ProcessCode::ABORT;

  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  auto isValidParticle = [&](const Data::SimParticle& p) {
    auto rho = std::hypot(p.position().x(), p.position().y());
    auto phi = std::atan2(p.momentum().y(), p.momentum().x());
    auto eta = std::atanh(p.momentum().z() / p.momentum().norm());
    auto pt  = std::hypot(p.momentum().x(), p.momentum().y());
    return within(rho, 0, m_cfg.rhoMax)
        and within(std::abs(p.position().z()), 0, m_cfg.absZMax)
        and within(phi, m_cfg.phiMin, m_cfg.phiMax)
        and within(eta, m_cfg.etaMin, m_cfg.etaMax)
        and within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax)
        and within(pt, m_cfg.ptMin, m_cfg.ptMax)
        and (m_cfg.keepNeutral or (p.q() != 0));
  };

  // copy selected vertices over to new collection
  size_t allParticles = 0;
  size_t selParticles = 0;

  for (const auto& vertex : *all) {

    allParticles += vertex.in.size();
    allParticles += vertex.out.size();

    Data::SimVertex<Data::SimParticle> sel;
    sel.position    = vertex.position;
    sel.timeStamp   = vertex.timeStamp;
    sel.processCode = vertex.processCode;

    // copy selected particles over
    std::copy_if(vertex.in.begin(),
                 vertex.in.end(),
                 std::back_inserter(sel.in),
                 isValidParticle);
    std::copy_if(vertex.out.begin(),
                 vertex.out.end(),
                 std::back_inserter(sel.out),
                 isValidParticle);

    // only retain vertex if it still contains particles
    if (not sel.in.empty() or not sel.out.empty()) {
      selParticles += sel.in.size();
      selParticles += sel.out.size();

      selected.push_back(std::move(sel));
    }
  }

  ACTS_DEBUG("event " << ctx.eventNumber << " selected " << selected.size()
                      << " from "
                      << all->size()
                      << " vertices");
  ACTS_DEBUG("event " << ctx.eventNumber << " selected " << selParticles
                      << " from "
                      << allParticles
                      << " particles");

  // write selected particles
  if (ctx.eventStore.add(m_cfg.output, std::move(selected))
      == ProcessCode::ABORT) {
    return ProcessCode::ABORT;
  }
  return ProcessCode::SUCCESS;
}
