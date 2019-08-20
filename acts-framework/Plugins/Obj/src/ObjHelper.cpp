// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Obj/ObjHelper.hpp"
#include <vector>

void
FW::Obj::writeVTN(std::ofstream&        stream,
                  VtnCounter&           vtnCounter,
                  double                scalor,
                  const Acts::Vector3D& vertex,
                  const std::string&    vtntype,
                  bool                  point)
{
  // in case you make a point
  unsigned int cp = 0;
  // the counter
  if (vtntype == "v") {
    ++vtnCounter.vcounter;
    cp = vtnCounter.vcounter;
  } else if (vtntype == "t") {
    ++vtnCounter.vtcounter;
    cp = vtnCounter.vtcounter;
  } else if (vtntype == "vn") {
    ++vtnCounter.ncounter;
    cp = vtnCounter.ncounter;
  } else
    return;

  // write out the vertex, texture vertex, normal
  stream << vtntype << " " << scalor * vertex.x() << " " << scalor * vertex.y()
         << " " << scalor * vertex.z() << '\n';
  // we create a point if needed
  if (point) stream << "p " << cp;
}

void
FW::Obj::constructVerticalFaces(std::ofstream&                   stream,
                                unsigned int                     start,
                                const std::vector<unsigned int>& vsides)
{
  // construct the vertical faces
  size_t       nsides = vsides.size();
  unsigned int sstart = start;
  for (auto vside : vsides) {
    if (vside) {
      // start streaming the side
      // all but the last
      if (start - sstart < nsides - 1) {
        stream << "f " << start << " " << start + 1 << " ";
        stream << start + nsides + 1 << " " << start + nsides;
      } else {
        stream << "f " << start << " " << sstart << " ";
        stream << sstart + nsides << " " << start + nsides;
      }
    }
    stream << '\n';
    // increase
    ++start;
  }
}

void
FW::Obj::writePlanarFace(std::ofstream&                     stream,
                         VtnCounter&                        vtnCounter,
                         double                             scalor,
                         const std::vector<Acts::Vector3D>& vertices,
                         double                             thickness,
                         const std::vector<unsigned int>&   vsides)
{
  // minimum 3 vertices needed
  if (vertices.size() < 3) return;
  // the first vertex
  unsigned int fvertex = vtnCounter.vcounter + 1;
  // lets create the normal vector first
  Acts::Vector3D sideOne = vertices[1] - vertices[0];
  Acts::Vector3D sideTwo = vertices[2] - vertices[1];
  Acts::Vector3D nvector(sideTwo.cross(sideOne).normalized());
  // write the normal vector
  writeVTN(stream, vtnCounter, scalor, nvector, "vn");
  // thickness or not thickness
  std::vector<int> sides     = {1};
  if (thickness != 0.) sides = {-1, 1};
  // now write all the vertices - this works w/wo thickness
  for (auto side : sides) {
    // save the current vertex counter
    unsigned int cvc = vtnCounter.vcounter;
    // loop over the sides
    for (auto v : vertices)
      writeVTN(stream,
               vtnCounter,
               scalor,
               v + (0.5 * side * thickness) * nvector,
               "v");
    // decide if you want to add texture
    std::string vtphr
        = "/";  // vtnCounter.vtcounter ? "/"+std::to_string(vtcounter) : "/";
    std::string ntphr = "/" + std::to_string(vtnCounter.ncounter);
    // now write the face
    stream << "f ";
    for (auto v : vertices) stream << ++cvc << vtphr << ntphr << " ";
    stream << '\n';
  }
  // now process the vertical sides
  constructVerticalFaces(stream, fvertex, vsides);
}

void
FW::Obj::writeTube(std::ofstream&           stream,
                   VtnCounter&              vtnCounter,
                   double                   scalor,
                   unsigned int             nSegments,
                   const Acts::Transform3D& transform,
                   double                   r,
                   double                   hZ,
                   double                   thickness)
{
  // flip along plus/minus and declare the faces
  std::vector<int> flip   = {-1, 1};
  std::vector<int> vfaces = {1, 2, 4, 3};
  // the number of phisteps
  double phistep = 2 * M_PI / nSegments;
  // make it twice if necessary
  std::vector<double> roffsets  = {0.};
  if (thickness != 0.) roffsets = {-0.5 * thickness, 0.5 * thickness};
  // now loop over the thickness and make an outer and inner
  unsigned int cvc   = vtnCounter.vcounter;
  size_t       iside = 0;
  for (auto t : roffsets) {
    size_t iphi = 0;
    // loop over phi steps
    for (; iphi < nSegments; ++iphi) {
      // currentPhi
      double phi = -M_PI + iphi * phistep;
      for (auto iflip : flip) {
        // create the vertex
        Acts::Vector3D point(transform * Acts::Vector3D((r + t) * cos(phi),
                                                        (r + t) * sin(phi),
                                                        iflip * hZ));
        // write the normal vector
        writeVTN(stream, vtnCounter, scalor, point, "v");
      }
    }
    // now create the faces
    iphi = 0;
    // side offset for faces
    unsigned int soff = 2 * iside * nSegments;
    for (; iphi < nSegments - 1; ++iphi) {
      // output to file
      stream << "f ";
      for (auto face : vfaces) stream << soff + cvc + (2 * iphi) + face << " ";
      stream << '\n';
    }
    // close the loop
    stream << "f " << soff + cvc + (2 * iphi) + 1 << " "
           << soff + cvc + (2 * iphi) + 2 << " " << soff + cvc + 2 << " "
           << soff + cvc + 1 << '\n';
    // new line at the end of the line
    stream << '\n';
    ++iside;
  }

  // construct the sides at the end when all vertices are done
  Acts::Vector3D nvectorSide = transform.rotation().col(2);
  // write the normal vector @todo flip sides
  writeVTN(stream, vtnCounter, scalor, nvectorSide, "vn");
  std::string ntphr = "//" + std::to_string(vtnCounter.ncounter);

  if (thickness != 0.) {
    // loop over the two sides
    for (iside = 0; iside < 2; ++iside) {
      // rest iphi
      size_t iphi = 0;
      for (; iphi < nSegments - 1; ++iphi) {
        stream << "f ";
        unsigned int base = cvc + (2 * iphi) + 1;
        stream << iside + base << ntphr << " ";
        stream << iside + base + 2 << ntphr << " ";
        stream << iside + base + (2 * nSegments) + 2 << ntphr << " ";
        stream << iside + base + (2 * nSegments) << ntphr << '\n';
      }
      // close the loop
      stream << "f ";
      stream << iside + cvc + (2 * iphi) + 1 << ntphr << " ";
      stream << iside + cvc + 1 << ntphr << " ";
      stream << iside + cvc + 1 + (2 * nSegments) << ntphr << " ";
      stream << iside + cvc + (2 * iphi) + 1 + (2 * nSegments) << ntphr << '\n';
    }
  }
}

// Bezier interpolation, see documentation
Acts::Vector3D
FW::Obj::calculateBezierPoint(double                t,
                              const Acts::Vector3D& p0,
                              const Acts::Vector3D& p1,
                              const Acts::Vector3D& p2,
                              const Acts::Vector3D& p3)
{
  double u   = 1. - t;
  double tt  = t * t;
  double uu  = u * u;
  double uuu = uu * u;
  double ttt = tt * t;

  Acts::Vector3D p = uuu * p0;  // first term
  p += 3 * uu * t * p1;         // second term
  p += 3 * u * tt * p2;         // third term
  p += ttt * p3;                // fourth term
  return p;
}
