// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Options/Pythia8Options.hpp"

#include "ACTFW/Generators/Pythia8ProcessGenerator.hpp"

void
FW::Options::addPythia8Options(boost::program_options::options_description& opt)
{
  using namespace boost::program_options;

  opt.add_options()("evg-cmsEnergy",
                    value<double>()->default_value(14000.),
                    "CMS value of the beam in [GeV].")(
      "evg-pdgBeam0",
      value<int>()->default_value(2212),
      "PDG number of beam 0 particles.")("evg-pdgBeam1",
                                         value<int>()->default_value(2212),
                                         "PDG number of beam 1 particles.")(
      "evg-hsProcess",
      value<std::string>()->default_value("HardQCD:all = on"),
      "The process string for the hard scatter event.")(
      "evg-puProcess",
      value<std::string>()->default_value("SoftQCD:all = on"),
      "The process string for the pile-up events.")(
      "evg-pileup",
      value<int>()->default_value(200),
      "Number of instantaneous pile-up events.")(
      "evg-vertex-xy-std",
      value<double>()->default_value(0.015),
      "Transverse vertex standard deviation in [mm].")(
      "evg-vertex-z-std",
      value<double>()->default_value(55.5),
      "Longitudinal vertex standard deviation in [mm].")(
      "evg-shuffle", bool_switch(), "Randomnly shuffle the vertex order.");
}

FW::EventGenerator::Config
FW::Options::readPythia8Options(const boost::program_options::variables_map& vm)
{
  Pythia8Generator::Config hardCfg;
  hardCfg.pdgBeam0  = vm["evg-pdgBeam0"].template as<int>();
  hardCfg.pdgBeam1  = vm["evg-pdgBeam1"].template as<int>();
  hardCfg.cmsEnergy = vm["evg-cmsEnergy"].template as<double>();
  hardCfg.settings  = {vm["evg-hsProcess"].template as<std::string>()};
  Pythia8Generator::Config pileupCfg;
  pileupCfg.pdgBeam0  = vm["evg-pdgBeam0"].template as<int>();
  pileupCfg.pdgBeam1  = vm["evg-pdgBeam1"].template as<int>();
  pileupCfg.cmsEnergy = vm["evg-cmsEnergy"].template as<double>();
  pileupCfg.settings  = {vm["evg-puProcess"].template as<std::string>()};

  auto vtxStdXY = vm["evg-vertex-xy-std"].template as<double>();
  auto vtxStdZ  = vm["evg-vertex-z-std"].template as<double>();

  EventGenerator::Config cfg;
  cfg.generators = {
      {FixedMultiplicityGenerator{1},
       GaussianVertexGenerator{vtxStdXY, vtxStdXY, vtxStdZ},
       Pythia8Generator::makeFunction(hardCfg)},
      {PoissonMultiplicityGenerator{
           static_cast<size_t>(vm["evg-pileup"].template as<int>())},
       GaussianVertexGenerator{vtxStdXY, vtxStdXY, vtxStdZ},
       Pythia8Generator::makeFunction(pileupCfg)},
  };
  cfg.shuffle = vm["evg-shuffle"].template as<bool>();

  return cfg;
}
