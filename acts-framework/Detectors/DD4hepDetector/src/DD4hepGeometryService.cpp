// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/DD4hepDetector/DD4hepGeometryService.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Plugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Tools/CylinderVolumeBuilder.hpp"
#include "Acts/Tools/CylinderVolumeHelper.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Tools/SurfaceArrayCreator.hpp"
#include "Acts/Tools/TrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/BinningType.hpp"

FW::DD4hep::DD4hepGeometryService::DD4hepGeometryService(
    const FW::DD4hep::DD4hepGeometryService::Config& cfg)
  : m_cfg(cfg), m_lcdd(), m_dd4hepGeometry(), m_trackingGeometry()
{
}

FW::DD4hep::DD4hepGeometryService::~DD4hepGeometryService()
{
  if (m_lcdd) m_lcdd->destroyInstance();
}

std::string
FW::DD4hep::DD4hepGeometryService::name() const
{
  return m_cfg.name;
}

FW::ProcessCode
FW::DD4hep::DD4hepGeometryService::buildDD4hepGeometry()
{
  m_lcdd = &(dd4hep::Detector::getInstance());
  for (auto& file : m_cfg.xmlFileNames) {
    m_lcdd->fromCompact(file.c_str());
  }
  m_lcdd->volumeManager();
  m_lcdd->apply("DD4hepVolumeManager", 0, 0);
  m_dd4hepGeometry = m_lcdd->world();

  return FW::ProcessCode::SUCCESS;
}

dd4hep::DetElement
FW::DD4hep::DD4hepGeometryService::dd4hepGeometry()
{
  if (!m_dd4hepGeometry) buildDD4hepGeometry();
  return m_dd4hepGeometry;
}

dd4hep::Detector*
FW::DD4hep::DD4hepGeometryService::DD4hepGeometryService::lcdd()
{
  if (!m_lcdd) buildDD4hepGeometry();
  return m_lcdd;
}

TGeoNode*
FW::DD4hep::DD4hepGeometryService::tgeoGeometry()
{
  if (!m_dd4hepGeometry) buildDD4hepGeometry();
  return m_dd4hepGeometry.placement().ptr();
}

FW::ProcessCode
FW::DD4hep::DD4hepGeometryService::buildTrackingGeometry()
{
  // set the tracking geometry
  m_trackingGeometry
      = std::move(Acts::convertDD4hepDetector(dd4hepGeometry(),
                                              m_cfg.lvl,
                                              m_cfg.bTypePhi,
                                              m_cfg.bTypeR,
                                              m_cfg.bTypeZ,
                                              m_cfg.envelopeR,
                                              m_cfg.envelopeZ,
                                              m_cfg.defaultLayerThickness,
                                              m_cfg.sortDetectors));
  return FW::ProcessCode::SUCCESS;
}

std::unique_ptr<const Acts::TrackingGeometry>
FW::DD4hep::DD4hepGeometryService::trackingGeometry()
{
  if (!m_trackingGeometry) buildTrackingGeometry();
  return std::move(m_trackingGeometry);
}
