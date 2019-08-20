// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/Pythia8ProcessGenerator.hpp"

#include <algorithm>
#include <random>

namespace {
struct FrameworkRndmEngine : public Pythia8::RndmEngine
{
  FW::RandomEngine& rng;

  FrameworkRndmEngine(FW::RandomEngine& rng_) : rng(rng_) {}
  double
  flat()
  {
    return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
  }
};
}  // namespace

std::function<std::vector<FW::Data::SimVertex<FW::Data::SimParticle>>(
    FW::RandomEngine&)>
FW::Pythia8Generator::makeFunction(const FW::Pythia8Generator::Config& cfg)
{
  auto gen = std::make_shared<Pythia8Generator>(cfg);
  return [=](RandomEngine& rng) { return (*gen)(rng); };
}

FW::Pythia8Generator::Pythia8Generator(const FW::Pythia8Generator::Config& cfg,
                                       Acts::Logging::Level level)
  : m_cfg(cfg)
  , m_logger(Acts::getDefaultLogger("Pythia8Generator", level))
  , m_pythia8("", false)
{
  // disable all output by default but allow reenable via config
  m_pythia8.settings.flag("Print:quiet", true);
  for (const auto& str : m_cfg.settings) {
    ACTS_VERBOSE("use Pythia8 setting '" << str << "'");
    m_pythia8.readString(str.c_str());
  }
  m_pythia8.settings.mode("Beams:idA", m_cfg.pdgBeam0);
  m_pythia8.settings.mode("Beams:idB", m_cfg.pdgBeam1);
  m_pythia8.settings.mode("Beams:frameType", 1);
  m_pythia8.settings.parm("Beams:eCM", m_cfg.cmsEnergy / Acts::units::_GeV);
  m_pythia8.init();
}

std::vector<FW::Data::SimVertex<FW::Data::SimParticle>>
FW::Pythia8Generator::operator()(FW::RandomEngine& rng)
{
  using namespace Data;

  // first vertex in list is the primary one at origin with time=0
  std::vector<SimVertex<SimParticle>> processes = {
      SimVertex<SimParticle>({0.0, 0.0, 0.0}),
  };

  // pythia8 is not thread safe and generation needs to be protected
  std::lock_guard<std::mutex> lock(m_pythia8Mutex);
  // use per-thread random engine also in pythia
  FrameworkRndmEngine rndmEngine(rng);
  m_pythia8.rndm.rndmEnginePtr(&rndmEngine);
  m_pythia8.next();

  // convert generated final state particles into acts format
  for (size_t ip = 0; ip < m_pythia8.event.size(); ++ip) {
    const auto& genParticle = m_pythia8.event[ip];

    // only interested in final, visible particles
    if (!genParticle.isFinal()) continue;
    if (!genParticle.isVisible()) continue;

    // define particle state
    // TODO check barcode and process code
    SimParticle particle(
        {
            genParticle.xProd() * Acts::units::_mm,
            genParticle.yProd() * Acts::units::_mm,
            genParticle.zProd() * Acts::units::_mm,
        },
        {
            genParticle.px() * Acts::units::_GeV,
            genParticle.py() * Acts::units::_GeV,
            genParticle.pz() * Acts::units::_GeV,
        },
        genParticle.m0(),
        genParticle.charge(),
        genParticle.id());

    if (genParticle.hasVertex()) {
      // either add to existing vertex w/ the same position or create new one
      // TODO can we do this w/p the manual search?
      auto it = std::find_if(processes.begin(),
                             processes.end(),
                             [&](const SimVertex<SimParticle>& vertex) {
                               return (vertex.position == particle.position());
                             });
      if (it == processes.end()) {
        processes.push_back(SimVertex<SimParticle>(particle.position()));
        processes.back().out.push_back(std::move(particle));
        ACTS_VERBOSE("created new secondary vertex "
                     << processes.back().position.transpose());
      } else {
        it->out.push_back(std::move(particle));
      }
    } else {
      // without a set vertex particles belong to the primary one
      processes.front().out.push_back(std::move(particle));
    }
  }
  return processes;
}
