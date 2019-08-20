// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/DD4hepDetector/DD4hepDetectorHelper.hpp"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Plugins/DD4hep/IActsExtension.hpp"
#include "DD4hep/DetFactoryHelper.h"

using namespace std;
using namespace dd4hep;

/**
 Constructor for a cylindrical barrel volume, possibly containing layers and the
 layers possibly containing modules. In Atlas style. This constructor was
 specifically introduced for the TrackML detector.
*/

static Ref_t
create_element(Detector& lcdd, xml_h xml, SensitiveDetector sens)
{
  xml_det_t x_det   = xml;
  string    detName = x_det.nameStr();
  // Make DetElement
  DetElement barrelDetector(detName, x_det.id());
  // add Extension to Detlement for the RecoGeometry
  Acts::ActsExtension::Config volConfig;
  volConfig.isBarrel             = true;
  Acts::ActsExtension* volumeExt = new Acts::ActsExtension(volConfig);
  barrelDetector.addExtension<Acts::IActsExtension>(volumeExt);
  // create Volume
  dd4hep::xml::Dimension x_detDim(x_det.dimensions());
  Volume                 barrelVolume(detName,
                      Tube(x_detDim.rmin(), x_detDim.rmax(), x_detDim.dz()),
                      lcdd.vacuum());  // vacuum at the moment change later
  barrelVolume.setVisAttributes(lcdd, x_detDim.visStr());
  // go trough possible layers
  size_t layerNumber = 0;

  for (xml_coll_t j(xml, _U(layer)); j; ++j) {
    xml_comp_t x_layer     = j;
    double     layerRmin   = x_layer.inner_r();
    double     layerRmax   = x_layer.outer_r();
    double     layerLength = x_layer.z();
    // Create Volume and DetElement for Layer
    string layerName = detName + _toString((int)layerNumber, "layer%d");
    Volume layerVolume(layerName,
                       Tube(layerRmin, layerRmax, layerLength),
                       lcdd.material(x_layer.materialStr()));
    DetElement layerDetector(barrelDetector, layerName, layerNumber);
    // Visualization
    layerVolume.setVisAttributes(lcdd, x_layer.visStr());
    // go trough possible modules
    if (x_layer.hasChild(_U(module))) {
      // the module describing the module dimensions
      xml_comp_t x_module         = x_layer.child(_U(module));
      double     modHalfLength    = x_module.length();
      double     modHalfWidth     = x_module.width();
      double     modHalfThickness = x_module.thickness();
      // informations for the placement and tilt in phi
      int    phiRepeat = x_module.repeat();
      double deltaphi  = 2. * M_PI / phiRepeat;
      double phiTilt   = x_module.phi_tilt();
      // informarions for he placement in z
      xml_comp_t x_slice  = x_layer.child(_U(slice));
      int        zRepeat  = x_slice.repeat();
      double     zOverlap = x_slice.dz();
      double     rOffset  = x_slice.offset();  // half offset

      // Create the module volume
      Volume modVolume("module",
                       Box(modHalfLength, modHalfWidth, modHalfThickness),
                       lcdd.material(x_module.materialStr()));

      // create the Acts::DigitizationModule (needed to do geometric
      // digitization) for all modules which have the same segmentation
      auto digiModule
          = FW::DD4hep::rectangleDigiModule(modHalfLength,
                                            modHalfWidth,
                                            modHalfThickness,
                                            sens.readout().segmentation());

      // Visualization
      modVolume.setVisAttributes(lcdd, x_module.visStr());
      size_t moduleNumber = 0;
      // Place the Modules in z
      // the distance between the modules in z
      double dz = (2 * modHalfLength - zOverlap);
      // the start value in z
      double startz = -((zRepeat - 1) * 0.5 * dz);
      // place in z
      for (int iz = 0; iz < zRepeat; iz++) {
        // to be added later to the module name
        string zname = _toString((int)iz, "z%d");
        // the radial position of the module
        double r = ((layerRmax + layerRmin) * 0.5);
        // alterning radial offset for each subsequent module
        if (iz % 2 == 0) {
          r += (0.5 * rOffset);
        } else {
          r -= (0.5 * rOffset);
        }
        // current z position
        double z = startz + iz * dz;
        // start phi position
        double minPhi = -M_PI + 0.5 * deltaphi / dd4hep::rad;
        // Place the modules in phi
        for (int iphi = 0; iphi < phiRepeat; ++iphi) {
          // the unique module name
          string module_name = zname + _toString((int)iphi, "module%d");
          // the phi position
          double phi = minPhi + deltaphi * iphi;
          // the position of the module within the layer
          Position trans(r * cos(phi), r * sin(phi), z);
          // Create the module DetElement
          DetElement moduleDetector(layerDetector, module_name, moduleNumber);
          // Set Sensitive Volmes sensitive
          if (x_module.isSensitive()) {
            modVolume.setSensitiveDetector(sens);
            // create and attach the extension with the shared digitzation
            // module
            Acts::ActsExtension* moduleExtension
                = new Acts::ActsExtension(digiModule);
            moduleDetector.addExtension<Acts::IActsExtension>(moduleExtension);
          }
          // Place Module Box Volumes in layer adding a tilt in phi
          PlacedVolume placedmodule = layerVolume.placeVolume(
              modVolume,
              Transform3D(RotationY(0.5 * M_PI) * RotationX(-phi - phiTilt),
                          trans));
          placedmodule.addPhysVolID("module", moduleNumber);
          // assign module DetElement to the placed module volume
          moduleDetector.setPlacement(placedmodule);
          ++moduleNumber;
        }
      }
    }
    // todo set granularity of layer material mapping and where material should
    // be mapped
    // hand over modules to ACTS
    Acts::ActsExtension::Config layConfig;
    layConfig.isLayer   = true;
    layConfig.envelopeR = 2. * Acts::units::_mm;
    layConfig.envelopeZ = 2. * Acts::units::_mm;  // maybe change later
    Acts::ActsExtension* detlayer = new Acts::ActsExtension(layConfig);
    layerDetector.addExtension<Acts::IActsExtension>(detlayer);
    // Place layer volume
    PlacedVolume placedLayer = barrelVolume.placeVolume(layerVolume);
    placedLayer.addPhysVolID("layer", layerNumber);
    // Assign layer DetElement to layer volume
    layerDetector.setPlacement(placedLayer);
    ++layerNumber;
  }
  // Place Volume
  Volume       mother_vol = lcdd.pickMotherVolume(barrelDetector);
  PlacedVolume placedTube = mother_vol.placeVolume(barrelVolume);
  placedTube.addPhysVolID("system", barrelDetector.id());
  barrelDetector.setPlacement(placedTube);

  return barrelDetector;
}

DECLARE_DETELEMENT(ACTS_GenericMLBarrel, create_element)
