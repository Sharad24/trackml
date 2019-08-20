// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/EventGenerator.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"

FW::EventGenerator::EventGenerator(const Config&        cfg,
                                   Acts::Logging::Level level)
  : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EventGenerator", level))
{
  if (m_cfg.output.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
  if (m_cfg.generators.empty()) {
    throw std::invalid_argument("No generators are configured");
  }
  if (!m_cfg.barcodeSvc) {
    throw std::invalid_argument("Missing barcode service");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers service");
  }
}

std::string
FW::EventGenerator::name() const
{
  return "EventGenerator";
}

size_t
FW::EventGenerator::numEvents() const
{
  return SIZE_MAX;
}

FW::ProcessCode
FW::EventGenerator::skip(size_t skip)
{
  // TODO 2018-03-13 msmk: how should skip be handled? ignore?
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::EventGenerator::read(AlgorithmContext ctx)
{
  std::vector<Data::SimVertex<Data::SimParticle>> event;

  auto   rng             = m_cfg.randomNumbers->spawnGenerator(ctx);
  size_t iGenerator      = 0;
  size_t iPrimary        = 0;  // primary vertex index within event
  size_t nTotalParticles = 0;  // total number of particles within event

  for (auto& generator : m_cfg.generators) {
    for (size_t n = generator.multiplicity(rng); 0 < n; --n) {

      // generate position and secondaries for this primary vertex
      auto vertex  = generator.vertex(rng);
      auto process = generator.process(rng);

      // modify secondaries to move everything to the primary vertex position
      size_t iSecondary = 0;  // secondary vertex index within current primary
      for (auto& secondaryVertex : process) {
        size_t iParticle = 0;  // particle index within secondary vertex

        // TODO use 4d vector in process directly
        Acts::Vector3D vertexPosition = vertex.head<3>();
        double         vertexTime     = vertex[3];
        secondaryVertex.position += vertexPosition;
        secondaryVertex.timeStamp += vertexTime;

        auto updateParticleInPlace = [&](Data::SimParticle& particle) {
          // generate new barcode and retain some existing information
          // TODO check if barcode components are correct
          auto generation  = m_cfg.barcodeSvc->generation(particle.barcode());
          auto processCode = m_cfg.barcodeSvc->process(particle.barcode());
          auto barcode     = m_cfg.barcodeSvc->generate(
              iPrimary, iParticle, generation, iSecondary, processCode);

          Acts::Vector3D position = particle.position() + vertexPosition;

          // TODO particle does not export timestamp?
          particle.place(position, barcode, vertexTime);
        };

        for (auto& particle : secondaryVertex.in) {
          updateParticleInPlace(particle);
          iParticle += 1;
        }
        for (auto& particle : secondaryVertex.out) {
          updateParticleInPlace(particle);
          iParticle += 1;
        }

        ACTS_VERBOSE("event " << ctx.eventNumber << " generator=" << iGenerator
                              << " primary="
                              << iPrimary
                              << " secondary="
                              << iSecondary
                              << " nparticles="
                              << iParticle);
        iSecondary += 1;
        nTotalParticles += iParticle;
      }

      // move all processes to the full event
      std::move(process.begin(), process.end(), std::back_inserter(event));

      iPrimary += 1;
    }
    iGenerator += 1;
  }
  if (m_cfg.shuffle) {
    std::shuffle(event.begin(), event.end(), rng);
  }

  ACTS_DEBUG("event " << ctx.eventNumber << " nprimaries=" << iPrimary
                      << " nparticles="
                      << nTotalParticles);

  // move generated event to the store
  if (ctx.eventStore.add(m_cfg.output, std::move(event))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }
  return FW::ProcessCode::SUCCESS;
}
