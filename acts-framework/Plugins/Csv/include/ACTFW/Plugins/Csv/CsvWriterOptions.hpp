// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include "ACTFW/Plugins/Csv/CsvSurfaceWriter.hpp"
#include "ACTFW/Plugins/Csv/CsvTrackingGeometryWriter.hpp"
#include "ACTFW/Utilities/Options.hpp"

namespace po = boost::program_options;

namespace FW {

namespace Options {

  // Common CSV writer option
  ///
  /// @tparam aopt_t Type of the options object (from BOOST)
  ///
  /// @param opt The options object, where string based options are attached
  template <typename aopt_t>
  void
  addCsvWriterOptions(aopt_t& opt)
  {
    opt.add_options()("csv-tg-fileheader",
                      po::value<std::string>()->default_value(""),
                      "The (optional) file header for the tracking geometry.")(
        "csv-tg-layerheader",
        po::value<std::string>()->default_value(""),
        "The (optional) header in front of layers.")(
        "csv-sf-fileheader",
        po::value<std::string>()->default_value(""),
        "The (optional) file header for the surface writer.")(
        "csv-sf-outputPrecission",
        po::value<int>()->default_value(6),
        "Floating number output precission.")(
        "csv-sf-outputScalor",
        po::value<double>()->default_value(1.),
        "Scale factor to be applied.")("csv-sf-outputBounds",
                                       po::value<bool>()->default_value(true),
                                       "Write the surface bounds to the file.")(
        "csv-sf-outputSensitive",
        po::value<bool>()->default_value(true),
        "Write sensitive surfaces.")("csv-sf-outputLayers",
                                     po::value<bool>()->default_value(true),
                                     "Write layer surfaces.");
  }

  /// Read the CSV options and return a Config file for the TrackingGeometry
  /// writing
  ///
  /// @tparam amap_t Type of the map object for reading out
  ///
  /// @param map Object for the options to be read out
  /// @param name Name of the TrackingGeometryWriter to be constructed
  template <typename amap_t>
  FW::Csv::CsvTrackingGeometryWriter::Config
  readCsvTrackingGeometryWriterConfig(const amap_t&      vm,
                                      const std::string& name
                                      = "CsvTrackingGeometryWriter")
  {
    FW::Csv::CsvTrackingGeometryWriter::Config objTgConfig(name,
                                                           Acts::Logging::INFO);
    objTgConfig.filePrefix = vm["csv-tg-fileheader"].template as<std::string>();
    objTgConfig.layerPrefix
        = vm["csv-tg-layerheader"].template as<std::string>();
    // Return the config object
    return objTgConfig;
  }

  /// Read the CSV options and return a Config file for the Surface
  /// writing
  ///
  /// @tparam amap_t Type of the map object for reading out
  ///
  /// @param map Object for the options to be read out
  /// @param name Name of the SurfaceWriter to be constructed
  template <typename amap_t>
  FW::Csv::CsvSurfaceWriter::Config
  readCsvSurfaceWriterConfig(const amap_t&      vm,
                             const std::string& name = "CSVSurfaceWriter")
  {
    FW::Csv::CsvSurfaceWriter::Config objSfConfig(name, Acts::Logging::INFO);

    objSfConfig.filePrefix = vm["csv-sf-fileheader"].template as<std::string>();
    objSfConfig.outputPrecision
        = vm["csv-sf-outputPrecission"].template as<int>();
    objSfConfig.outputScalor = vm["csv-sf-outputScalor"].template as<double>();
    objSfConfig.outputBounds = vm["csv-sf-outputBounds"].template as<bool>();
    objSfConfig.outputSensitive
        = vm["csv-sf-outputSensitive"].template as<bool>();
    objSfConfig.outputLayerSurface
        = vm["csv-sf-outputLayers"].template as<bool>();
    // Return the config object
    return objSfConfig;
  }

}  // namespace Options
}  // namespace FW
