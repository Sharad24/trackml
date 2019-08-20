// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//-------------------------------------------------------------------------------------
// 3 passive layers
//-------------------------------------------------------------------------------------
auto material = Acts::Material(100., 300., 9.012, 4., 1.848e-3);
Acts::PassiveLayerBuilder::Config sphenixConfig;
sphenixConfig.layerIdentification = "SPhenix";
sphenixConfig.centralLayerRadii   = {23.4, 31.5, 39., 60.0, 80.0, 100.0, 120.0};
sphenixConfig.centralLayerHalflengthZ
    = {270., 270., 270., 500., 500., 500., 500.};
sphenixConfig.centralLayerThickness = {1., 1., 1., 1., 1., 1., 1.};
sphenixConfig.centralLayerMaterial
    = {material, material, material, material, material, material, material};

double rmin   = 300.;
double rmax   = 780.;
int    layers = 78;

double rstep = (rmax - rmin) / layers;
for (int il = 0; il < layers; ++il) {
  double rcurrent = rmin + il * rstep;
  sphenixConfig.centralLayerRadii.push_back(rcurrent);
  sphenixConfig.centralLayerHalflengthZ.push_back(1055.);
  sphenixConfig.centralLayerThickness.push_back(0.01);
  sphenixConfig.centralLayerMaterial.push_back(material);
}

auto sphenixBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
    sphenixConfig,
    Acts::getDefaultLogger("SphenixLayerBuilder", layerLLevel));
// create the volume for the beam pipe
Acts::CylinderVolumeBuilder::Config sphenixvolConfig;
sphenixvolConfig.trackingVolumeHelper = cylinderVolumeHelper;
sphenixvolConfig.volumeName           = "Sphenix";
sphenixvolConfig.buildToRadiusZero    = false;
sphenixvolConfig.layerBuilder         = sphenixBuilder;
sphenixvolConfig.volumeSignature      = 0;
auto sphenixVolumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
    sphenixvolConfig,
    Acts::getDefaultLogger("SPhenixLayerBuilder", volumeLLevel));

// add to the detector builds
volumeBuilders.push_back(sphenixVolumeBuilder);