// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include <fstream>
#include <ios>
#include <stdexcept>
#include "ACTFW/Barcode/Barcode.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"

FW::Csv::CsvParticleWriter::CsvParticleWriter(
    const FW::Csv::CsvParticleWriter::Config& cfg,
    Acts::Logging::Level                      level)
  : ParticleWriter(cfg.collection, "CsvParticleWriter", level), m_cfg(cfg)
{
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
}

FW::ProcessCode
FW::Csv::CsvParticleWriter::writeT(
    const FW::AlgorithmContext&           ctx,
    const std::vector<Data::SimVertex<>>& vertices)
{
  std::string pathOs = perEventFilepath(
      m_cfg.outputDir, m_cfg.outputFileName, ctx.eventNumber);
  std::ofstream os(pathOs, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + pathOs + "' to write");
  }

  // do we have the hits per particle written out ?
  bool hppPresent = m_cfg.hitsPerParticleCollection != "";

  const std::map<barcode_type, size_t>* hitsPerParticle = nullptr;
  if (hppPresent
      && ctx.eventStore.get(m_cfg.hitsPerParticleCollection, hitsPerParticle)
          == ProcessCode::ABORT) {
    throw std::ios_base::failure(
        "Could not retrieve hits/particle reference map.");
  }

  // write csv header
  os << "particle_id,";
  os << "vx,vy,vz,";
  os << "px,py,pz,";
  os << "q";
  if (hppPresent) os << ",nhits";
  os << '\n';

  // write one line per particle
  os << std::setprecision(m_cfg.outputPrecision);
  for (auto& vertex : vertices) {
    for (auto& particle : vertex.outgoing()) {
      os << particle.barcode() << ",";
      os << particle.position().x() << ",";
      os << particle.position().y() << ",";
      os << particle.position().z() << ",";
      os << particle.momentum().x() << ",";
      os << particle.momentum().y() << ",";
      os << particle.momentum().z() << ",";
      os << particle.q();
      // add the hits per particle
      if (hitsPerParticle) {
        auto hppEntry = hitsPerParticle->find(particle.barcode());
        int  nhits    = 0;
        if (hppEntry != hitsPerParticle->end()) {
          nhits = hppEntry->second;
          os << "," << nhits;
        }
      }
      os << '\n';
    }
  }

  return ProcessCode::SUCCESS;
}
