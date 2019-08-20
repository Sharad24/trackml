// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include <iostream>
#include "ACTFW/Writers/IWriterT.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Surfaces/Surface.hpp"

FW::Obj::ObjTrackingGeometryWriter::ObjTrackingGeometryWriter(
    const FW::Obj::ObjTrackingGeometryWriter::Config& cfg)
  : FW::IWriterT<Acts::TrackingGeometry>(), m_cfg(cfg)
{
}

std::string
FW::Obj::ObjTrackingGeometryWriter::name() const
{
  return m_cfg.name;
}

FW::ProcessCode
FW::Obj::ObjTrackingGeometryWriter::write(
    const Acts::TrackingGeometry& tGeometry)
{
  ACTS_DEBUG(">>Obj: Writer for TrackingGeometry object called.");
  // get the world volume
  auto world = tGeometry.highestTrackingVolume();
  if (world) write(*world);
  // return the success code
  return FW::ProcessCode::SUCCESS;
}

/// process this volume
void
FW::Obj::ObjTrackingGeometryWriter::write(const Acts::TrackingVolume& tVolume)
{
  ACTS_DEBUG(">>Obj: Writer for TrackingVolume object called.");
  // get the confined layers and process them
  if (tVolume.confinedLayers()) {
    ACTS_VERBOSE(">>Obj: Layers are present, process them.");
    // loop over the layers
    for (auto layer : tVolume.confinedLayers()->arrayObjects()) {
      // we jump navigation layers
      if (layer->layerType() == Acts::navigation) continue;
      // get the volume name
      const std::string& volumeName = tVolume.volumeName();
      // find the right surfacewriter
      std::shared_ptr<ObjSurfaceWriter> surfaceWriter = nullptr;
      for (auto writer : m_cfg.surfaceWriters) {
        // get name and writer
        auto writerName = writer->name();
        // and break
        ACTS_VERBOSE(">>Obj: The writer name is: " << writerName);
        ACTS_VERBOSE(">>Obj: The volume name is: " << volumeName);
        if (volumeName.find(writerName) != std::string::npos) {
          // asign the writer
          surfaceWriter = writer;
          // break the loop
          break;
        }
      }
      // bail out if you have no surface writer
      if (!surfaceWriter) return;
      // layer prefix
      surfaceWriter->write(m_cfg.layerPrefix);
      // try to write the material surface as well
      if (layer->surfaceRepresentation().associatedMaterial())
        surfaceWriter->write(layer->surfaceRepresentation());
      // the the approaching surfaces and check if they have material
      if (layer->approachDescriptor()) {
        // loop over the contained Surfaces
        for (auto& cSurface : layer->approachDescriptor()->containedSurfaces())
          if (cSurface->associatedMaterial()) surfaceWriter->write(*cSurface);
      }
      // check for sensitive surfaces
      if (layer->surfaceArray() && surfaceWriter) {
        ACTS_VERBOSE(">>Obj: There are "
                     << layer->surfaceArray()->surfaces().size()
                     << " surfaces.");
        // surfaces
        // surfaceWriter->write(m_cfg.sensitiveGroupPrefix);
        // loop over the surface
        for (auto& surface : layer->surfaceArray()->surfaces()) {
          if (surface
              && (surfaceWriter->write(*surface)) == FW::ProcessCode::ABORT)
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
