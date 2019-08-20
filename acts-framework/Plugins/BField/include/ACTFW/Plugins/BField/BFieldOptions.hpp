// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_OPTIONS_BFIELDOPTIONS_HPP
#define ACTFW_OPTIONS_BFIELDOPTIONS_HPP

#include <iostream>
#include <utility>
#include "ACTFW/Plugins/BField/BFieldUtils.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  // common bfield options, with a bf prefix
  template <class AOPT>
  void
  addBFieldOptions(AOPT& opt)
  {
    opt.add_options()(
        "bf-map",
        po::value<std::string>()->default_value(""),
        "Set this string to point to the bfield source file."
        "That can either be a '.txt', a '.csv' or a '.root' file. "
        "Omit for a constant magnetic field.")(
        "bf-name",
        po::value<std::string>()->default_value("bField"),
        "In case your field map file is given in root format, please specify "
        "the "
        "name of the TTree.")(
        "bf-gridpoints",
        po::value<size_t>()->default_value(100000),
        "Estimate of number of grid points, "
        "needed for allocation, only for txt and csv files.")(
        "bf-lscalor",
        po::value<double>()->default_value(1.),
        "The default unit for the grid "
        "points is mm. In case the grid points of your field map has another "
        "unit, please set  the scalor to mm.")(
        "bs-bscalor",
        po::value<double>()->default_value(1.),
        "The default unit for the magnetic field values is Tesla. In case the "
        "grid points of your field map has another unit, please set  the "
        "scalor "
        "to [T].")(
        "bf-rz",
        po::value<bool>()->default_value(false),
        "Please set this flag to true, if your grid points and your "
        "magnetic field values are given in 'rz'. The default is 'xyz'.")(
        "bf-foctant",
        po::value<bool>()->default_value(false),
        "Please set this flag to true, if your field map is only given for the "
        "first octant/quadrant and should be symmetrically created for all "
        "other "
        "octants/quadrants.")(
        "bf-values",
        po::value<read_range>()->multitoken()->default_value({0., 0., 0.}),
        "In case no magnetic field map is handed over. A constant magnetic "
        "field will be created automatically. The values can be set with this "
        "options. Please hand over the coordinates in cartesian coordinates: "
        "{Bx,By,Bz} in Tesla.");
  }

  // create the bfield maps
  template <class AMAP>
  std::pair<std::shared_ptr<Acts::InterpolatedBFieldMap>,
            std::shared_ptr<Acts::ConstantBField>>
  readBField(const AMAP& vm)
  {
    std::string bfieldmap = "constfield";

    enum BFieldMapType { constant = 0, root = 1, text = 2 };

    int bfieldmaptype = constant;
    if (vm.count("bf-map") && vm["bf-map"].template as<std::string>() != "") {
      bfieldmap = vm["bf-map"].template as<std::string>();
      std::cout << "- read in magnetic field map: "
                << vm["bf-map"].template as<std::string>() << std::endl;
      if (bfieldmap.find(".root") != std::string::npos) {
        std::cout << "- root format for magnetic field detected" << std::endl;
        bfieldmaptype = root;
      } else if (bfieldmap.find(".txt") != std::string::npos
                 || bfieldmap.find(".csv") != std::string::npos) {
        std::cout << "- txt format for magnetic field detected" << std::endl;
        bfieldmaptype = text;
      } else {
        std::cout << "- magnetic field format could not be detected";
        std::cout << " use '.root', '.txt', or '.csv'." << std::endl;
        return std::pair<std::shared_ptr<Acts::InterpolatedBFieldMap>,
                         std::shared_ptr<Acts::ConstantBField>>(nullptr,
                                                                nullptr);
      }
    }
    if (bfieldmaptype == text && vm.count("bf-gridpoints")) {
      std::cout << "- number of points set to: "
                << vm["bf-gridpoints"].template as<size_t>() << std::endl;
    }
    double lscalor = 1.;
    if (bfieldmaptype != constant && vm.count("bf-lscalor")) {
      lscalor = vm["bf-lscalor"].template as<double>();
      std::cout << "- length scalor to mm set to: " << lscalor << std::endl;
    }
    double bscalor = 1.;
    if (vm.count("bf-bscalor")) {
      bscalor = vm["bf-bscalor"].template as<double>();
      std::cout << "- BField (scalor to/in) Tesla set to: " << bscalor
                << std::endl;
    }
    if (bfieldmaptype != constant && vm["bf-rz"].template as<bool>())
      std::cout << "- BField map is given in 'rz' coordiantes." << std::endl;
    else if (bfieldmaptype != constant)
      std::cout << "- BField map is given in 'xyz' coordiantes." << std::endl;

    if (bfieldmaptype != constant && vm["bf-foctant"].template as<bool>()) {
      std::cout
          << "- Only the first octant/quadrant is given, bField map will be "
             "symmetrically created for all other octants/quadrants"
          << std::endl;
    }

    // Declare the mapper
    Acts::concept::AnyFieldLookup<> mapper;
    double                          lengthUnit = lscalor * Acts::units::_mm;
    double                          BFieldUnit = bscalor * Acts::units::_T;

    // set the mapper - foort
    if (bfieldmaptype == root) {
      if (vm["bf-rz"].template as<bool>()) {
        mapper = FW::BField::root::fieldMapperRZ(
            [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
              return (binsRZ.at(1) * nBinsRZ.at(0) + binsRZ.at(0));
            },
            vm["bf-map"].template as<std::string>(),
            vm["bf-name"].template as<std::string>(),
            lengthUnit,
            BFieldUnit,
            vm["bf-foctant"].template as<bool>());
      } else {
        mapper = FW::BField::root::fieldMapperXYZ(
            [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
              return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2))
                      + binsXYZ.at(1) * nBinsXYZ.at(2)
                      + binsXYZ.at(2));
            },
            vm["bf-map"].template as<std::string>(),
            vm["bf-name"].template as<std::string>(),
            lengthUnit,
            BFieldUnit,
            vm["bf-foctant"].template as<bool>());
      }
    } else if (bfieldmaptype == text) {
      if (vm["bf-rz"].template as<bool>()) {
        mapper = FW::BField::txt::fieldMapperRZ(
            [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
              return (binsRZ.at(1) * nBinsRZ.at(0) + binsRZ.at(0));
            },
            vm["bf-map"].template as<std::string>(),
            lengthUnit,
            BFieldUnit,
            vm["bf-gridpoints"].template as<size_t>(),
            vm["bf-foctant"].template as<bool>());
      } else {
        mapper = FW::BField::txt::fieldMapperXYZ(
            [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
              return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2))
                      + binsXYZ.at(1) * nBinsXYZ.at(2)
                      + binsXYZ.at(2));
            },
            vm["bf-map"].template as<std::string>(),
            lengthUnit,
            BFieldUnit,
            vm["bf-gridpoints"].template as<size_t>(),
            vm["bf-foctant"].template as<bool>());
      }
    }
    Acts::InterpolatedBFieldMap::Config config;
    config.scale  = 1.;
    config.mapper = std::move(mapper);

    std::shared_ptr<Acts::InterpolatedBFieldMap> bField
        = bfieldmaptype != constant
        ? std::make_shared<Acts::InterpolatedBFieldMap>(std::move(config))
        : nullptr;

    // No bfield map is handed over
    // get the constant bField values
    auto bFieldValues = vm["bf-values"].template as<read_range>();
    if (bFieldValues.size() != 3) {
      throw std::invalid_argument(
          "- The values handed over for the constant magnetic field "
          "have wrong dimension. Needs to have 3 dimension. Please "
          "hand over the coordinates in cartesian coordinates: "
          "{Bx,By,Bz} in Tesla.");
    }
    // Create the constant magnetic field
    std::shared_ptr<Acts::ConstantBField> cField
        = std::make_shared<Acts::ConstantBField>(
            bFieldValues.at(0) * Acts::units::_T,
            bFieldValues.at(1) * Acts::units::_T,
            bFieldValues.at(2) * Acts::units::_T);

    return std::pair<std::shared_ptr<Acts::InterpolatedBFieldMap>,
                     std::shared_ptr<Acts::ConstantBField>>(bField, cField);
  }
}
}

#endif  // ACTFW_OPTIONS_BFIELDOPTIONS_HPP
