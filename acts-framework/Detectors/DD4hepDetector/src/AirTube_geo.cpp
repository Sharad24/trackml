// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This DD4hep Constructor is strongly inspired by (adds the ActsExtension)
// https://github.com/AIDASoft/DD4hep/tree/master/examples/SimpleDetector

// $Id: $
//====================================================================
//  // Simple tube filled with air
//  // used for tracking purposes ...
//
//--------------------------------------------------------------------
//
//  Author     : F.Gaede
//
//====================================================================
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace dd4hep;

static Ref_t
create_element(Detector& lcdd, xml_h e, SensitiveDetector /* sens */)
{
  xml_det_t   x_det = e;
  std::string name  = x_det.nameStr();

  DetElement airTube(name, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension::Config volConfig;
  volConfig.isBeampipe           = true;
  Acts::ActsExtension* detvolume = new Acts::ActsExtension(volConfig);
  airTube.addExtension<Acts::IActsExtension>(detvolume);

  PlacedVolume pv;

  // ----- read xml ----------------------

  xml_dim_t dim = x_det.dimensions();

  double inner_r    = dim.rmin();
  double outer_r    = dim.rmax();
  double z_half     = dim.zhalf();
  double tube_thick = outer_r - inner_r;

  //--------------------------------------

  Tube tubeSolid(inner_r, outer_r, z_half);

  Volume tube_vol(
      name + "_inner_cylinder_air", tubeSolid, lcdd.material("Air"));

  /*  Vector3D ocyl(inner_r + 0.5 * tube_thick, 0., 0.);

    VolCylinder cylSurf(tube_vol,
                        SurfaceType(SurfaceType::Helper),
                        0.5 * tube_thick,
                        0.5 * tube_thick,
                        ocyl);

    volSurfaceList(airTube)->push_back(cylSurf);*/

  //--------------------------------------

  Volume mother = lcdd.pickMotherVolume(airTube);

  pv = mother.placeVolume(tube_vol);

  pv.addPhysVolID("system", x_det.id());

  airTube.setPlacement(pv);

  return airTube;
}

DECLARE_DETELEMENT(AirTube, create_element)
