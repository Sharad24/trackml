// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Empty.hpp"

#include <stdexcept>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/GeometryID.hpp>

#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FW::EmptyReconstructionAlgorithm::EmptyReconstructionAlgorithm(
    const FW::EmptyReconstructionAlgorithm::Config& cfg,
    Acts::Logging::Level                            logLevel)
  : BareAlgorithm("EmptyReconstructionAlgorithm", logLevel), m_cfg(cfg)
{
  if (m_cfg.spacePointCollection.empty()) {
    throw std::invalid_argument("Missing input space points collection");
  }
}

FW::ProcessCode
FW::EmptyReconstructionAlgorithm::execute(FW::AlgorithmContext ctx) const
{
  ACTS_INFO("empty reconstruction on event " << ctx.eventNumber);

  // initialize point struct
  struct point {
    double x, y, z;
    point() {}
    point(double x, double y, double z) : x(x),y(y),z(z) {}
    inline point operator-(const point&p) {
      return point(x-p.x, y-p.y, z-p.z);
    }
    inline point operator+(const point&p) {
      return point(x+p.x, y+p.y, z+p.z);
    }
    inline double operator*(const point&p) {
      return x*p.x+y*p.y+z*p.z;
    }
    inline point operator*(double f) {
      return point(x*f, y*f, z*f);
    }
  };

  // initialize vector of points to store hits
  std::vector<point> hits;

  // initialize character array and file name
  char file[1000];
  sprintf(file, "/tmp/event000000000-hits.csv");

  FILE*fp = fopen(file, "r");
  if (!fp) {
    ACTS_INFO("couldn't open hits\n");
    exit(1);
  }
  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);

  hits.push_back(point(0,0,0));

  ACTS_INFO("Reading data")
  while (1) {
    // char out_str[50];
    long long hit_id;
    double tx, ty, tz;
    int volume_id, layer_id, module_id;
    if (fscanf(fp, "%lld,%lf,%lf,%lf,%d,%d,%d", &hit_id, &tx, &ty, &tz, &volume_id, &layer_id, &module_id) == -1) break;
    hits.push_back(point(tx, ty, tz));
    // sprintf(out_str, "%lld, %lf, %lf, %lf", hit_id, tx, ty, tz);
    /// ACTS_INFO(out_str);
  }
  ACTS_INFO("Read hits data");
  fclose(fp);

  // initialize structure for cells data
  // pair<pair<ch0, ch1>, value>
  std::vector<std::pair<std::pair<int, int>, double> > hit_cells[200000];

  //initialize filename and file pointers
  char file2[1000];
  sprintf(file2, "/tmp/event000000000-details.csv");

  FILE*fp2 = fopen(file2, "r");
  if (!fp2) {
    printf("couldn't open cells\n");
    exit(1);
  }
  ACTS_INFO("Reading cells data");
  char tmpstr2[1000];
  int tmp2 = fscanf(fp2, "%s", tmpstr);
  int hit_id, ch0, ch1;
  double value;
  while (fscanf(fp2, "%d,%d,%d,%lf", &hit_id, &ch0, &ch1, &value) == 4) {
    hit_cells[hit_id].push_back(std::make_pair(std::make_pair(ch0, ch1), value));
  }
  fclose(fp2);
  ACTS_INFO("Read cells data");
  // Read the input space points
  const DetectorData<geo_id_value, Acts::Vector3D>* spacePoints;
  if (ctx.eventStore.get(m_cfg.spacePointCollection, spacePoints)
      == FW::ProcessCode::ABORT) {
    ACTS_WARNING("missing space point input");
    return FW::ProcessCode::ABORT;
  }
  ACTS_INFO("Loaded input points for event " << ctx.eventNumber);
  for (auto& volumeData : *spacePoints) {
    for (auto& layerData : volumeData.second) {
      for (auto& moduleData : layerData.second) {
        ACTS_INFO("volume " << volumeData.first << " layer " << layerData.first
                            << " module " << moduleData.first << " "
                            << moduleData.second.size() << " space points");
        for (auto& cluster : moduleData.second) {
          // TODO do something w/ the space points and write out track
          // candidates NOTE space points are grouped by the geometry id
        }
      }
    }
  }

  return ProcessCode::SUCCESS;
}
