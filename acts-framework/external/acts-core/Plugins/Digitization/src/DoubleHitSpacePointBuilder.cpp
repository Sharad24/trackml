// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Digitization/DoubleHitSpacePointBuilder.hpp"
#include <cmath>
#include <limits>
#include "Acts/Utilities/Helpers.hpp"

///
/// @note Used abbreviation: "Strip Detector Element" -> SDE
///
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::SpacePointBuilder(
    Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::
        DoubleHitSpacePointConfig cfg)
  : m_cfg(std::move(cfg))
{
}

double
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::differenceOfClustersChecked(
    const Acts::Vector3D& pos1,
    const Acts::Vector3D& pos2) const
{
  // Check if measurements are close enough to each other
  if ((pos1 - pos2).norm() > m_cfg.diffDist) {
    return -1.;
  }

  // Calculate the angles of the vectors
  double phi1, theta1, phi2, theta2;
  phi1   = VectorHelpers::phi(pos1 - m_cfg.vertex);
  theta1 = VectorHelpers::theta(pos1 - m_cfg.vertex);
  phi2   = VectorHelpers::phi(pos2 - m_cfg.vertex);
  theta2 = VectorHelpers::theta(pos2 - m_cfg.vertex);

  // Calculate the squared difference between the theta angles
  double diffTheta2 = (theta1 - theta2) * (theta1 - theta2);
  if (diffTheta2 > m_cfg.diffTheta2) {
    return -1.;
  }

  // Calculate the squared difference between the phi angles
  double diffPhi2 = (phi1 - phi2) * (phi1 - phi2);
  if (diffPhi2 > m_cfg.diffPhi2) {
    return -1.;
  }

  // Return the squared distance between both vector
  return diffTheta2 + diffPhi2;
}

Acts::Vector2D
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::localCoords(
    const Acts::PlanarModuleCluster& cluster) const
{
  // Local position information
  auto           par = cluster.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

Acts::Vector3D
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::globalCoords(
    const Acts::PlanarModuleCluster& cluster) const
{
  // Receive corresponding surface
  auto& clusterSurface = cluster.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(localCoords(cluster), mom, pos);

  return pos;
}

void
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::makeClusterPairs(
    const std::vector<const PlanarModuleCluster*>& clustersFront,
    const std::vector<const PlanarModuleCluster*>& clustersBack,
    std::vector<std::pair<const Acts::PlanarModuleCluster*,
                          const Acts::PlanarModuleCluster*>>& clusterPairs)
    const
{
  // Return if no clusters are given in a vector
  if (clustersFront.empty() || clustersBack.empty()) {
    return;
  }

  // Declare helper variables
  double       currentDiff;
  double       diffMin;
  unsigned int clusterMinDist;

  // Walk through all clusters on both surfaces
  for (unsigned int iClustersFront = 0; iClustersFront < clustersFront.size();
       iClustersFront++) {
    // Set the closest distance to the maximum of double
    diffMin = std::numeric_limits<double>::max();
    // Set the corresponding index to an element not in the list of clusters
    clusterMinDist = clustersBack.size();
    for (unsigned int iClustersBack = 0; iClustersBack < clustersBack.size();
         iClustersBack++) {
      // Calculate the distances between the hits
      currentDiff = differenceOfClustersChecked(
          globalCoords(*(clustersFront[iClustersFront])),
          globalCoords(*(clustersBack[iClustersBack])));
      // Store the closest clusters (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin        = currentDiff;
        clusterMinDist = iClustersBack;
      }
    }

    // Store the best (=closest) result
    if (clusterMinDist < clustersBack.size()) {
      std::pair<const Acts::PlanarModuleCluster*,
                const Acts::PlanarModuleCluster*>
          clusterPair;
      clusterPair = std::make_pair(clustersFront[iClustersFront],
                                   clustersBack[clusterMinDist]);
      clusterPairs.push_back(clusterPair);
    }
  }
}

std::pair<Acts::Vector3D, Acts::Vector3D>
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::endsOfStrip(
    const Acts::PlanarModuleCluster& cluster) const
{
  // Calculate the local coordinates of the cluster
  const Acts::Vector2D local = localCoords(cluster);

  // Receive the binning
  auto segment = dynamic_cast<const Acts::CartesianSegmentation*>(
      &(cluster.digitizationModule()->segmentation()));
  auto& binData     = segment->binUtility().binningData();
  auto& boundariesX = binData[0].boundaries();
  auto& boundariesY = binData[1].boundaries();

  // Search the x-/y-bin of the cluster
  size_t binX = binData[0].searchLocal(local);
  size_t binY = binData[1].searchLocal(local);

  Acts::Vector2D topLocal, bottomLocal;

  if (boundariesX[binX + 1] - boundariesX[binX]
      < boundariesY[binY + 1] - boundariesY[binY]) {
    // Set the top and bottom end of the strip in local coordinates
    topLocal = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
                boundariesY[binY + 1]};
    bottomLocal
        = {(boundariesX[binX] + boundariesX[binX + 1]) / 2, boundariesY[binY]};
  } else {
    // Set the top and bottom end of the strip in local coordinates
    topLocal
        = {boundariesX[binX], (boundariesY[binY] + boundariesY[binY + 1]) / 2};
    bottomLocal = {boundariesX[binX + 1],
                   (boundariesY[binY] + boundariesY[binY + 1]) / 2};
  }

  // Calculate the global coordinates of the top and bottom end of the strip
  Acts::Vector3D topGlobal, bottomGlobal, mom;  // mom is a dummy variable
  const auto*    sur = &cluster.referenceSurface();
  sur->localToGlobal(topLocal, mom, topGlobal);
  sur->localToGlobal(bottomLocal, mom, bottomGlobal);

  // Return the top and bottom end of the strip in global coordinates
  return std::make_pair(topGlobal, bottomGlobal);
}

double
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::calcPerpProj(
    const Acts::Vector3D& a,
    const Acts::Vector3D& c,
    const Acts::Vector3D& q,
    const Acts::Vector3D& r) const
{
  /// This approach assumes that no vertex is available. This option aims to
  /// approximate the space points from cosmic data.
  /// The underlying assumption is that the best point is given by the closest
  /// distance between both lines describing the SDEs.
  /// The point x on the first SDE is parametrized as a + lambda0 * q with the
  /// top end a of the strip and the vector q = a - b(ottom end of the strip).
  /// An analogous parametrization is performed of the second SDE with y = c +
  /// lambda1 * r.
  /// x get resolved by resolving lambda0 from the condition that |x-y| is the
  /// shortest distance between two skew lines.
  Acts::Vector3D ac    = c - a;
  double         qr    = q.dot(r);
  double         denom = q.dot(q) - qr * qr;

  // Check for numerical stability
  if (fabs(denom) > 10e-7) {
    // Return lambda0
    return (ac.dot(r) * qr - ac.dot(q) * r.dot(r)) / denom;
  }
  // lambda0 is in the interval [-1,0]. This return serves as error check.
  return 1.;
}

bool
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::recoverSpacePoint(
    Acts::SpacePointBuilder<DoubleHitSpacePoint>::SpacePointParameters& spaPoPa)
    const
{
  /// Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  if (m_cfg.stripLengthGapTolerance <= 0.) {
    return false;
  }
  spaPoPa.qmag = spaPoPa.q.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spaPoPa.limitExtended
      = spaPoPa.limit + m_cfg.stripLengthGapTolerance / spaPoPa.qmag;
  // Check if m is just slightly outside
  if (fabs(spaPoPa.m) > spaPoPa.limitExtended) {
    return false;
  }
  // Calculate n if not performed previously
  if (spaPoPa.n == 0.) {
    spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs);
  }
  // Check if n is just slightly outside
  if (fabs(spaPoPa.n) > spaPoPa.limitExtended) {
    return false;
  }

  /// The following code considers an overshoot of m and n in the same direction
  /// of their SDE. The term "overshoot" represents the amount of m or n outside
  /// its regular interval (-1, 1).
  /// It calculates which overshoot is worse. In order to compare both, the
  /// overshoot in n is projected onto the first surface by considering the
  /// normalized projection of r onto q.
  /// This allows a rescaling of the overshoot. The worse overshoot will be set
  /// to +/-1, the parameter with less overshoot will be moved towards 0 by the
  /// worse overshoot.
  /// In order to treat both SDEs equally, the rescaling eventually needs to be
  /// performed several times. If these shifts allows m and n to be in the
  /// limits, the space point can be stored.
  /// @note This shift can be understood as a shift of the particle's
  /// trajectory. This is leads to a shift of the vertex. Since these two points
  /// are treated independently from other measurement, it is also possible to
  /// consider this as a change in the slope of the particle's trajectory. The
  /// would also move the vertex position.

  // Calculate the scaling factor to project lengths of the second SDE on the
  // first SDE
  double secOnFirstScale
      = spaPoPa.q.dot(spaPoPa.r) / (spaPoPa.qmag * spaPoPa.qmag);
  // Check if both overshoots are in the same direction
  if (spaPoPa.m > 1. && spaPoPa.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spaPoPa.m - 1.;
    double nOvershoot
        = (spaPoPa.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spaPoPa.m -= biggerOvershoot;
    spaPoPa.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // Check if both overshoots are in the same direction
  if (spaPoPa.m < -1. && spaPoPa.n < -1.) {
    // Calculate the overshoots
    double mOvershoot = -(spaPoPa.m + 1.);
    double nOvershoot
        = -(spaPoPa.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spaPoPa.m += biggerOvershoot;
    spaPoPa.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // No solution could be found
  return false;
}

void
Acts::SpacePointBuilder<Acts::DoubleHitSpacePoint>::calculateSpacePoints(
    const std::vector<std::pair<const Acts::PlanarModuleCluster*,
                                const Acts::PlanarModuleCluster*>>&
                                            clusterPairs,
    std::vector<Acts::DoubleHitSpacePoint>& spacePoints) const
{

  /// Source of algorithm: Athena, SiSpacePointMakerTool::makeSCT_SpacePoint()

  Acts::SpacePointBuilder<DoubleHitSpacePoint>::SpacePointParameters spaPoPa;

  // Walk over every found candidate pair
  for (const auto& cp : clusterPairs) {

    // Calculate the ends of the SDEs
    const auto& ends1 = endsOfStrip(*(cp.first));
    const auto& ends2 = endsOfStrip(*(cp.second));

    /// The following algorithm is meant for finding the position on the first
    /// strip if there is a corresponding cluster on the second strip. The
    /// resulting point is a point x on the first surfaces. This point is
    /// along a line between the points a (top end of the strip)
    /// and b (bottom end of the strip). The location can be parametrized as
    /// 	2 * x = (1 + m) a + (1 - m) b
    /// as function of the scalar m. m is a parameter in the interval
    /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
    /// from the vertex to the cluster on the second strip y is needed to be a
    /// multiple k of the vector from vertex to the hit on the first strip x.
    /// As a consequence of this demand y = k * x needs to be on the
    /// connecting line between the top (c) and bottom (d) end of
    /// the second strip. If both clusters correspond to each other, the
    /// condition
    /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
    /// needs to be fulfilled. Inserting the first equation into this
    /// equation leads to the condition for m as given in the following
    /// algorithm and therefore to the calculation of x.
    /// The same calculation can be repeated for y. Its corresponding
    /// parameter will be named n.

    spaPoPa.q = ends1.first - ends1.second;
    spaPoPa.r = ends2.first - ends2.second;

    // Fast skipping if a perpendicular projection should be used
    double resultPerpProj;
    if (m_cfg.usePerpProj) {
      resultPerpProj
          = calcPerpProj(ends1.first, ends2.first, spaPoPa.q, spaPoPa.r);
      if (resultPerpProj <= 0.) {
        Acts::DoubleHitSpacePoint sp;
        sp.clusterFront = cp.first;
        sp.clusterBack  = cp.second;
        sp.spacePoint   = ends1.first + resultPerpProj * spaPoPa.q;
        spacePoints.push_back(std::move(sp));
        continue;
      }
    }

    spaPoPa.s  = ends1.first + ends1.second - 2 * m_cfg.vertex;
    spaPoPa.t  = ends2.first + ends2.second - 2 * m_cfg.vertex;
    spaPoPa.qs = spaPoPa.q.cross(spaPoPa.s);
    spaPoPa.rt = spaPoPa.r.cross(spaPoPa.t);
    spaPoPa.m  = -spaPoPa.s.dot(spaPoPa.rt) / spaPoPa.q.dot(spaPoPa.rt);

    // Set the limit for the parameter
    if (spaPoPa.limit == 1. && m_cfg.stripLengthTolerance != 0.) {
      spaPoPa.limit = 1. + m_cfg.stripLengthTolerance;
    }

    // Check if m and n can be resolved in the interval (-1, 1)
    if (fabs(spaPoPa.m) <= spaPoPa.limit
        && fabs(spaPoPa.n
                = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs))
            <= spaPoPa.limit) {
      // Store the space point
      Acts::DoubleHitSpacePoint sp;
      sp.clusterFront = cp.first;
      sp.clusterBack  = cp.second;
      sp.spacePoint
          = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
      spacePoints.push_back(std::move(sp));
    } else {
      /// If this point is reached then it was not possible to resolve both
      /// points such that they are on their SDEs
      /// The following code treats a possible recovery of points resolved
      /// slightly outside of the SDE.
      /// @note This procedure is an indirect variation of the vertex
      /// position.
      // Check if a recovery the point(s) and store them if successful
      if (recoverSpacePoint(spaPoPa)) {
        Acts::DoubleHitSpacePoint sp;
        sp.clusterFront = cp.first;
        sp.clusterBack  = cp.second;
        sp.spacePoint
            = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);
        spacePoints.push_back(std::move(sp));
      }
    }
  }
}
