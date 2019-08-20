// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Csv/CsvSurfaceWriter.hpp"
#include <ios>
#include <iostream>
#include <stdexcept>
#include "Acts/Layers/Layer.hpp"
#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/GeometryID.hpp"

FW::Csv::CsvSurfaceWriter::CsvSurfaceWriter(
    const FW::Csv::CsvSurfaceWriter::Config& cfg)
  : FW::IWriterT<Acts::Surface>(), m_cfg(cfg)
{
  // Validate the configuration
  if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  } else if (m_cfg.name.empty()) {
    throw std::invalid_argument("Missing algorithm name");
  } else if (!m_cfg.outputStream) {
    throw std::invalid_argument("Missing output stream");
  }

  // Write down the file prefix
  if (m_cfg.filePrefix != "") (*m_cfg.outputStream) << m_cfg.filePrefix << '\n';
  (*m_cfg.outputStream) << "volume_id,layer_id,module_id,";
  (*m_cfg.outputStream) << "cx,cy,cz,";
  (*m_cfg.outputStream) << "rot_xu,rot_xv,rot_xw,";
  (*m_cfg.outputStream) << "rot_yu,rot_yv,rot_yw,";
  (*m_cfg.outputStream) << "rot_zu,rot_zv,rot_zw";
  if (m_cfg.outputBounds) {
    (*m_cfg.outputStream)
        << ",module_t,module_minhx,module_maxhx,module_hy,pitchX,pitchY";
  }
  (*m_cfg.outputStream) << '\n';
  (*m_cfg.outputStream) << std::setprecision(m_cfg.outputPrecision);
}

std::string
FW::Csv::CsvSurfaceWriter::name() const
{
  return m_cfg.name;
}

FW::ProcessCode
FW::Csv::CsvSurfaceWriter::write(const Acts::Surface& surface)
{
  std::lock_guard<std::mutex> lock(m_write_mutex);

  // check
  ACTS_DEBUG(">>Csv: Writer for Surface object called.");

  // let's get the bounds & the transform
  const Acts::SurfaceBounds& surfaceBounds = surface.bounds();
  auto                       sCenter       = surface.center();
  auto                       sTransform    = surface.transform();

  // Get the layer geo id information
  auto geoID       = surface.geoID();
  auto volumeID    = geoID.value(Acts::GeometryID::volume_mask);
  auto layerID     = geoID.value(Acts::GeometryID::layer_mask);
  auto sensitiveID = geoID.value(Acts::GeometryID::sensitive_mask);

  // write configurations
  if ((sensitiveID && m_cfg.outputSensitive)
      || (!sensitiveID && m_cfg.outputLayerSurface)) {

    // surface indentification
    (*m_cfg.outputStream) << volumeID << "," << layerID << "," << sensitiveID
                          << ",";
    (*m_cfg.outputStream) << sCenter.x() << "," << sCenter.y() << ","
                          << sCenter.z() << ",";
    (*m_cfg.outputStream) << sTransform(0, 0) << "," << sTransform(0, 1) << ","
                          << sTransform(0, 2) << ",";
    (*m_cfg.outputStream) << sTransform(1, 0) << "," << sTransform(1, 1) << ","
                          << sTransform(1, 2) << ",";
    (*m_cfg.outputStream) << sTransform(2, 0) << "," << sTransform(2, 1) << ","
                          << sTransform(2, 2);

    // dynamic_cast to PlanarBounds
    const Acts::PlanarBounds* planarBounds
        = dynamic_cast<const Acts::PlanarBounds*>(&surfaceBounds);
    // only continue if the cast worked
    if (planarBounds && m_cfg.outputSensitive && m_cfg.outputBounds) {
      ACTS_VERBOSE(">>Csv: Writing out a PlaneSurface ");
      // get thickness and value store
      if (surface.associatedDetectorElement()) {
        auto detElement = dynamic_cast<const Acts::IdentifiedDetectorElement*>(
            surface.associatedDetectorElement());
        // get thickness and bounds
        double thickness = detElement->thickness();
        (*m_cfg.outputStream) << "," << thickness << ",";
        // get the values from the bound value store
        auto bValues = surface.bounds().valueStore();
        if (bValues.size() == 2) {
          (*m_cfg.outputStream) << bValues[0] << ",";
          (*m_cfg.outputStream) << bValues[0] << ",";
          (*m_cfg.outputStream) << bValues[1] << ",";
        } else if (bValues.size() == 3) {
          for (auto& bv : bValues) (*m_cfg.outputStream) << bv << ",";
        }
        // get the ditigization module
        if (detElement->digitizationModule()) {
          auto dModule = detElement->digitizationModule();
          // dynamic_cast to CartesianSegmentation
          const Acts::CartesianSegmentation* cSegmentation
              = dynamic_cast<const Acts::CartesianSegmentation*>(
                  &(dModule->segmentation()));
          if (cSegmentation) {
            auto pitch = cSegmentation->pitch();
            (*m_cfg.outputStream) << pitch.first << ",";
            (*m_cfg.outputStream) << pitch.second;
          }
        }
      }
    }
    // break the line
    (*m_cfg.outputStream) << '\n';

    // check if you have layer and check what your have
    // dynamic cast to CylinderBounds work the same
    const Acts::CylinderBounds* cylinderBounds
        = dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
    if (cylinderBounds && m_cfg.outputLayerSurface) {
      ACTS_VERBOSE(">>Csv: Writing out a CylinderSurface with r = "
                   << cylinderBounds->r());
      (*m_cfg.outputStream) << '\n';
    }

    ////dynamic cast to RadialBounds or disc bounds work the same
    const Acts::RadialBounds* radialBounds
        = dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
    if (radialBounds && m_cfg.outputLayerSurface) {
      ACTS_VERBOSE(">>Csv: Writing out a DiskSurface at z = "
                   << sTransform.translation().z());
      // name the object
      // we use the tube writer in the other direction
      // double rMin      = radialBounds->rMin();
      // double rMax      = radialBounds->rMax();
      // double thickness = rMax - rMin;
      (*m_cfg.outputStream) << '\n';
    }
  }
  // return success
  return FW::ProcessCode::SUCCESS;
}
