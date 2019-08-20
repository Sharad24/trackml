// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Digitization/DigitizationAlgorithm.hpp"

#include <iostream>
#include <stdexcept>

#include "ACTFW/Barcode/Barcode.hpp"
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "Acts/Detector/DetectorElementBase.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleStepper.hpp"
#include "Acts/Plugins/Digitization/Segmentation.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

FW::DigitizationAlgorithm::DigitizationAlgorithm(
    const FW::DigitizationAlgorithm::Config& cfg,
    Acts::Logging::Level                     level)
  : FW::BareAlgorithm("DigitizationAlgorithm", level), m_cfg(cfg)
{
  if (m_cfg.simulatedHitCollection.empty()) {
    throw std::invalid_argument("Missing input hits collection");
  } else if (m_cfg.spacePointCollection.empty()) {
    throw std::invalid_argument("Missing output space points collection");
  } else if (m_cfg.clusterCollection.empty()) {
    throw std::invalid_argument("Missing output clusters collection");
  } else if (!m_cfg.randomNumberSvc) {
    throw std::invalid_argument("Missing random numbers service");
  } else if (!m_cfg.planarModuleStepper) {
    throw std::invalid_argument("Missing planar module stepper");
  }
}

FW::ProcessCode
FW::DigitizationAlgorithm::execute(FW::AlgorithmContext context) const
{
  // Prepare the input data collection
  const FW::DetectorData<geo_id_value, Data::SimHit<Data::SimParticle>>* simHits
      = nullptr;

  // Read it from the event store
  if (context.eventStore.get(m_cfg.simulatedHitCollection, simHits)
      == FW::ProcessCode::ABORT)
    return FW::ProcessCode::ABORT;

  ACTS_DEBUG("Retrieved hit data '" << m_cfg.simulatedHitCollection
                                    << "' from event store.");

  // Prepare the output data: Clusters
  FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster> planarClusters;

  // Prepare the second output data : Truth SpacePoints (for debugging)
  FW::DetectorData<geo_id_value, Acts::Vector3D> spacePoints;

  // now digitise
  for (auto& vData : (*simHits)) {
    auto volumeKey = vData.first;
    ACTS_DEBUG("- Processing Volume Data collection for volume with ID "
               << volumeKey);
    for (auto& lData : vData.second) {
      auto layerKey = lData.first;
      ACTS_DEBUG("-- Processing Layer Data collection for layer with ID "
                 << layerKey);
      for (auto& sData : lData.second) {
        auto moduleKey = sData.first;
        ACTS_DEBUG("-- Processing Module Data collection for module with ID "
                   << moduleKey);
        // get the hit parameters
        for (auto& hit : sData.second) {
          // get the surface
          const Acts::Surface& hitSurface = (*hit.surface);
          // get the associated particle from the hit
          const Data::SimParticle*              hitParticle = &hit.particle;
          std::vector<const Data::SimParticle*> hitParticles{hitParticle};
          // get the DetectorElement
          auto hitDetElement
              = dynamic_cast<const Acts::IdentifiedDetectorElement*>(
                  hitSurface.associatedDetectorElement());
          if (hitDetElement) {
            // get the digitization module
            auto hitDigitizationModule = hitDetElement->digitizationModule();
            if (hitDigitizationModule) {
              // get the lorentz angle
              double lorentzAngle = hitDigitizationModule->lorentzAngle();
              double thickness    = hitDetElement->thickness();
              double lorentzShift = thickness * tan(lorentzAngle);
              lorentzShift *= -(hitDigitizationModule->readoutDirection());
              // parameters
              auto invTransfrom = hitSurface.transform().inverse();
              // local intersection / direction
              Acts::Vector3D localIntersect3D(invTransfrom * hit.position);
              Acts::Vector2D localIntersection(localIntersect3D.x(),
                                               localIntersect3D.y());
              Acts::Vector3D localDirection(invTransfrom.linear()
                                            * hit.direction);
              // now calculate the steps through the silicon
              std::vector<Acts::DigitizationStep> dSteps
                  = m_cfg.planarModuleStepper->cellSteps(*hitDigitizationModule,
                                                         localIntersection,
                                                         localDirection);
              // everything under threshold or edge effects
              if (!dSteps.size()) continue;
              /// let' create a cluster - centroid method
              double localX    = 0.;
              double localY    = 0.;
              double totalPath = 0.;
              // the cells to be used
              std::vector<Acts::DigitizationCell> usedCells;
              usedCells.reserve(dSteps.size());
              // loop over the steps
              for (auto dStep : dSteps) {
                // @todo implement smearing
                localX += dStep.stepLength * dStep.stepCellCenter.x();
                localY += dStep.stepLength * dStep.stepCellCenter.y();
                totalPath += dStep.stepLength;
                usedCells.push_back(
                    Acts::DigitizationCell(dStep.stepCell.channel0,
                                           dStep.stepCell.channel1,
                                           dStep.stepLength));
              }
              // divide by the total path
              localX /= totalPath;
              localX += lorentzShift;
              localY /= totalPath;

              // get the segmentation & find the corresponding cell id
              const Acts::Segmentation& segmentation
                  = hitDigitizationModule->segmentation();
              auto           binUtility = segmentation.binUtility();
              Acts::Vector2D localPosition(localX, localY);
              // @todo remove unneccesary conversion
              size_t bin0          = binUtility.bin(localPosition, 0);
              size_t bin1          = binUtility.bin(localPosition, 1);
              size_t binSerialized = binUtility.serialize({{bin0, bin1, 0}});

              // the covariance is currently set to 0.
              Acts::ActsSymMatrixD<2> cov;
              cov << 0.05, 0., 0.05, 0.;

              // create the geometry based Idnetifier
              Acts::GeometryID geoID(0);
              geoID.add(volumeKey, Acts::GeometryID::volume_mask);
              geoID.add(layerKey, Acts::GeometryID::layer_mask);
              geoID.add(moduleKey, Acts::GeometryID::sensitive_mask);
              geoID.add(binSerialized, Acts::GeometryID::channel_mask);

              // create the planar cluster
              Acts::PlanarModuleCluster pCluster(
                  hitSurface.getSharedPtr(),
                  Identifier(Identifier::identifier_type(geoID.value()),
                             hitParticles),
                  std::move(cov),
                  localX,
                  localY,
                  std::move(usedCells));

              // insert into the cluster map
              FW::Data::insert(planarClusters,
                               volumeKey,
                               layerKey,
                               moduleKey,
                               std::move(pCluster));

            }  // hit moulde proection
          }    // hit element protection
        }      // hit loop
      }        // moudle loop
    }          // layer loop
  }            // volume loop

  // write the SpacePoints to the EventStore
  if (context.eventStore.add(m_cfg.spacePointCollection, std::move(spacePoints))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }
  // write the clusters to the EventStore
  if (context.eventStore.add(m_cfg.clusterCollection, std::move(planarClusters))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }

  return FW::ProcessCode::SUCCESS;
}
