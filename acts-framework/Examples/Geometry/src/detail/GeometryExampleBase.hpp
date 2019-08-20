// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>
#include <string>
#include <vector>
#include "ACTFW/Common/CommonOptions.hpp"
#include "ACTFW/Common/GeometryOptions.hpp"
#include "ACTFW/Common/OutputOptions.hpp"
#include "ACTFW/Plugins/Csv/CsvSurfaceWriter.hpp"
#include "ACTFW/Plugins/Csv/CsvTrackingGeometryWriter.hpp"
#include "ACTFW/Plugins/Csv/CsvWriterOptions.hpp"
#include "ACTFW/Plugins/Obj/ObjSurfaceWriter.hpp"
#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "ACTFW/Plugins/Obj/ObjWriterOptions.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"

namespace po = boost::program_options;

template <typename geometry_options_t, typename geometry_getter_t>
int
processGeometry(int                argc,
                char*              argv[],
                geometry_options_t geometryOptions,
                geometry_getter_t  trackingGeometry)
{

  // Declare the supported program options.
  po::options_description desc("Allowed options");
  // Add the standard/common options
  FW::Options::addCommonOptions<po::options_description>(desc);
  // Add options for the Obj writing
  FW::Options::addObjWriterOptions<po::options_description>(desc);
  // Add the geometry options
  FW::Options::addGeometryOptions<po::options_description>(desc);
  // Add the output options
  FW::Options::addOutputOptions<po::options_description>(desc);
  // Add specific options for this geometry
  geometryOptions(desc);

  po::variables_map vm;
  // Get all options from contain line and store it into the map
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  // Print help if requested
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  // Now read the standard options
  auto logLevel  = FW::Options::readLogLevel<po::variables_map>(vm);
  auto tGeometry = trackingGeometry(vm);

  // The detectors
  read_strings subDetectors = vm["geo-subdetectors"].as<read_strings>();

  Acts::Logging::Level surfaceLogLevel
      = Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  Acts::Logging::Level layerLogLevel
      = Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  Acts::Logging::Level volumeLogLevel
      = Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  // ---------------------------------------------------------------------------------
  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();

  // OBJ output
  if (vm["output-obj"].as<bool>()) {
    // The writers
    std::vector<std::shared_ptr<FW::Obj::ObjSurfaceWriter>> subWriters;
    std::vector<std::shared_ptr<std::ofstream>>             subStreams;
    // Loop and create the obj output writers per defined sub detector
    for (auto sdet : subDetectors) {
      // Sub detector stream
      auto        sdStream = std::shared_ptr<std::ofstream>(new std::ofstream);
      std::string sdOutputName = FW::joinPaths(outputDir, sdet + ".obj");
      sdStream->open(sdOutputName);
      // Object surface writers
      FW::Obj::ObjSurfaceWriter::Config sdObjWriterConfig
          = FW::Options::readObjSurfaceWriterConfig(vm, sdet, surfaceLogLevel);
      sdObjWriterConfig.outputStream = sdStream;
      auto sdObjWriter
          = std::make_shared<FW::Obj::ObjSurfaceWriter>(sdObjWriterConfig);
      // Push back
      subWriters.push_back(sdObjWriter);
      subStreams.push_back(sdStream);
    }

    // Configure the tracking geometry writer
    auto tgObjWriterConfig = FW::Options::readObjTrackingGeometryWriterConfig(
        vm, "ObjTrackingGeometryWriter", volumeLogLevel);

    tgObjWriterConfig.surfaceWriters = subWriters;
    auto tgObjWriter = std::make_shared<FW::Obj::ObjTrackingGeometryWriter>(
        tgObjWriterConfig);

    // Write the tracking geometry object
    tgObjWriter->write(*tGeometry);

    // Close the output streams
    for (auto sStreams : subStreams) {
      sStreams->close();
    }
  }

  // CSV output
  if (vm["output-csv"].as<bool>()) {

    auto        csvStream = std::shared_ptr<std::ofstream>(new std::ofstream);
    std::string csvOutputName = "Detector.csv";
    csvStream->open(csvOutputName);

    FW::Csv::CsvSurfaceWriter::Config sfCsvWriterConfig
        = FW::Options::readCsvSurfaceWriterConfig(vm, "CsvSurfaceWriter");
    sfCsvWriterConfig.outputStream = csvStream;
    auto sfCsvWriter
        = std::make_shared<FW::Csv::CsvSurfaceWriter>(sfCsvWriterConfig);

    // Configure the tracking geometry writer
    FW::Csv::CsvTrackingGeometryWriter::Config tgCsvWriterConfig
        = FW::Options::readCsvTrackingGeometryWriterConfig(
            vm, "CsvTrackingGeometryWriter");
    tgCsvWriterConfig.surfaceWriter = sfCsvWriter;
    tgCsvWriterConfig.layerPrefix   = "\n";
    // The tracking geometry writer
    auto tgCsvWriter = std::make_shared<FW::Csv::CsvTrackingGeometryWriter>(
        tgCsvWriterConfig);

    // Write the tracking geometry object
    tgCsvWriter->write(*(tGeometry.get()));

    // Close the file
    csvStream->close();
  }

  return 0;
}
