// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <fstream>
#include <ios>
#include <stdexcept>

#include <Acts/Plugins/Digitization/PlanarModuleCluster.hpp>
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/EventData/SimIdentifier.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Plugins/Csv/CsvPlanarClusterWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"

FW::Csv::CsvPlanarClusterWriter::CsvPlanarClusterWriter(
    const FW::Csv::CsvPlanarClusterWriter::Config& cfg,
    Acts::Logging::Level                           level)
  : Base(cfg.collection, "CsvPlanarClusterWriter", level), m_cfg(cfg)
{
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  }
}

FW::ProcessCode
FW::Csv::CsvPlanarClusterWriter::writeT(
    const AlgorithmContext& ctx,
    const FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster>& clusters)
{
  // open per-event hits file
  std::string pathHits
      = perEventFilepath(m_cfg.outputDir, "hits.csv", ctx.eventNumber);
  std::ofstream osHits(pathHits, std::ofstream::out | std::ofstream::trunc);
  if (!osHits) {
    throw std::ios_base::failure("Could not open '" + pathHits + "' to write");
  }

  // open per-event details file for the hit details
  std::string pathDetails
      = perEventFilepath(m_cfg.outputDir, "details.csv", ctx.eventNumber);
  std::ofstream osDetails(pathDetails,
                          std::ofstream::out | std::ofstream::trunc);
  if (!osDetails) {
    throw std::ios_base::failure("Could not open '" + pathHits + "' to write");
  }

  // open per-event truth file
  std::string pathTruth
      = perEventFilepath(m_cfg.outputDir, "truth.csv", ctx.eventNumber);
  std::ofstream osTruth(pathTruth, std::ofstream::out | std::ofstream::trunc);
  if (!osTruth) {
    throw std::ios_base::failure("Could not open '" + pathTruth + "' to write");
  }

  size_t skipped_hits = 0;

  // write csv hits header
  osHits << "hit_id,";
  osHits << "x,y,z,";
  osHits << "volume_id,layer_id,module_id" << '\n';
  osHits << std::setprecision(m_cfg.outputPrecision);

  // write csv hit detials header
  osDetails << "hit_id,ch0,ch1,value" << '\n';
  osDetails << std::setprecision(m_cfg.outputPrecision);

  // write csv truth headers
  osTruth << "hit_id,";
  osTruth << "particle_id,";
  osTruth << "tx,ty,tz,";
  osTruth << "tpx,tpy,tpz\n";

  size_t hitId = 0;
  for (auto& volumeData : clusters) {
    for (auto& layerData : volumeData.second) {
      for (auto& moduleData : layerData.second) {
        for (auto& cluster : moduleData.second) {
          hitId += 1;
          // local cluster information
          auto           parameters = cluster.parameters();
          Acts::Vector2D local(parameters[Acts::ParDef::eLOC_0],
                               parameters[Acts::ParDef::eLOC_1]);
          Acts::Vector3D pos(0, 0, 0);
          Acts::Vector3D mom(1, 1, 1);
          // transform local into global position information
          cluster.referenceSurface().localToGlobal(local, mom, pos);

          // write hit information
          osHits << hitId << ",";
          osHits << pos.x() << "," << pos.y() << "," << pos.z() << ",";
          osHits << volumeData.first << ",";
          osHits << layerData.first << ",";
          osHits << moduleData.first << '\n';

          // append cell information
          auto cells = cluster.digitizationCells();
          for (auto& cell : cells) {
            osDetails << hitId << "," << cell.channel0 << "," << cell.channel1
                      << "," << cell.data << '\n';
          }
          /// Hit identifier
          auto hitIdentifier = cluster.sourceLink();
          // write hit-particle truth association
          // each hit can have multiple particles, e.g. in a dense environment
          for (auto& sPartilce : hitIdentifier.truthParticles()) {
            // positon
            const Acts::Vector3D& sPosition = sPartilce->position();
            const Acts::Vector3D& sMomentum = sPartilce->position();
            osTruth << hitId << "," << sPartilce->barcode() << ",";
            osTruth << sPosition.x() << "," << sPosition.y() << ","
                    << sPosition.z() << ',';
            osTruth << sMomentum.x() << "," << sMomentum.y() << ","
                    << sMomentum.z() << '\n';
          }
        }
      }
    }
  }

  ACTS_VERBOSE(
      "Number of skipped hits from being written out : " << skipped_hits);

  return FW::ProcessCode::SUCCESS;
}
