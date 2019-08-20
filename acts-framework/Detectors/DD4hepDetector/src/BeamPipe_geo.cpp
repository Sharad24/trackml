// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

static Ref_t
create_element(Detector& lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det    = e;
  string    det_name = x_det.nameStr();
  // Make DetElement
  DetElement beamtube(det_name, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension::Config volConfig;
  volConfig.isBeampipe           = true;
  Acts::ActsExtension* detvolume = new Acts::ActsExtension(volConfig);
  beamtube.addExtension<Acts::IActsExtension>(detvolume);

  dd4hep::xml::Dimension x_det_dim(x_det.dimensions());
  Tube   tube_shape(x_det_dim.rmin(), x_det_dim.rmax(), x_det_dim.z());
  Volume tube_vol(
      det_name, tube_shape, lcdd.air());  // air at the moment change later
  tube_vol.setVisAttributes(lcdd, x_det_dim.visStr());
  // Place Volume
  Volume       mother_vol = lcdd.pickMotherVolume(beamtube);
  PlacedVolume placedTube = mother_vol.placeVolume(tube_vol);
  placedTube.addPhysVolID("tube", beamtube.id());
  beamtube.setPlacement(placedTube);

  return beamtube;
}

DECLARE_DETELEMENT(BeamPipe, create_element)
