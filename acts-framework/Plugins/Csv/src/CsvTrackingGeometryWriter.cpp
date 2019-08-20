// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Csv/CsvTrackingGeometryWriter.hpp"
#include <iostream>
#include <sstream>
#include "ACTFW/Writers/IWriterT.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"

FW::Csv::CsvTrackingGeometryWriter::CsvTrackingGeometryWriter(
    const FW::Csv::CsvTrackingGeometryWriter::Config& cfg)
  : FW::IWriterT<Acts::TrackingGeometry>(), m_cfg(cfg)
{
}

std::string
FW::Csv::CsvTrackingGeometryWriter::name() const
{
  return m_cfg.name;
}

FW::ProcessCode
FW::Csv::CsvTrackingGeometryWriter::write(
    const Acts::TrackingGeometry& tGeometry)
{
  ACTS_DEBUG(">>Csv: Writer for TrackingGeometry object called.");
  // get the world volume
  auto world = tGeometry.highestTrackingVolume();
  if (world) write(*world);
  // return the success code
  return FW::ProcessCode::SUCCESS;
}

/// process this volume
void
FW::Csv::CsvTrackingGeometryWriter::write(const Acts::TrackingVolume& tVolume)
{
  ACTS_DEBUG(">>Csv: Writer for TrackingVolume object called.");
  // get the confined layers and process them
  if (tVolume.confinedLayers()) {
    ACTS_VERBOSE(">>Csv: Layers are present, process them.");
    // loop over the layers
    for (auto layer : tVolume.confinedLayers()->arrayObjects()) {
      // we jump navigation layers
      if (layer->layerType() == Acts::navigation) continue;
      // find the right surfacewriter
      auto surfaceWriter = m_cfg.surfaceWriter;
      // bail out if you have no surface writer
      if (!surfaceWriter) return;
      // layer prefix
      surfaceWriter->write(m_cfg.layerPrefix);
      // check for sensitive surfaces
      if (layer->surfaceArray() && surfaceWriter) {
        // the current module thickness
        std::vector<double> cValues;
        // loop over the surface
        for (auto surface : layer->surfaceArray()->surfaces()) {
          if (surface
              && surfaceWriter->write(*surface) == FW::ProcessCode::ABORT)
            return;
        }
      }
    }
  }
  // get the confined volumes and step down the hierarchy
  if (tVolume.confinedVolumes()) {
    // loop over the volumes and write what they have
    for (auto volume : tVolume.confinedVolumes()->arrayObjects()) {
      write(*volume.get());
    }
  }
}
