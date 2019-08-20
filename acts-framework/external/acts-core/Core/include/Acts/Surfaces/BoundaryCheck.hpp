// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundaryCheck.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <cfloat>
#include <cmath>
#include <iterator>
#include <vector>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @class BoundaryCheck
///
/// The BoundaryCheck class provides boundary checks and distance calculations
/// for aligned box-like and polygonal boundaries on local surfaces.
/// Different types of boundary checks are supported and are transparently
/// selected when calling the `isInside(...)` and `distance(...)` methods:
///
/// -   Hard checks w/o any tolerances
/// -   Tolerance-based checks in one or in both local coordinates
/// -   Chi2-based checks based on a covariance matrix. Non-vanishing
///     correlations are correctly taken into account.
///
/// With a defined covariance matrix, the closest point and the distance are
/// not defined along the usual Euclidean metric, but by the Mahalanobis
/// distance induced by the the covariance.
class BoundaryCheck
{
public:
  /// Construct either hard cut in both dimensions or no cut at all.
  BoundaryCheck(bool check);

  /// Construct a tolerance based check.
  ///
  /// @param checkLocal0 Boolean directive to check coordinate 0
  /// @param checkLocal1 Boolean directive to check coordinate 1
  /// @param tolerance0 Tolerance along coordinate 0
  /// @param tolerance1 Tolerance along coordinate 1
  BoundaryCheck(bool   checkLocal0,
                bool   checkLocal1,
                double tolerance0 = 0,
                double tolerance1 = 0);

  /// Construct a chi2-based check.
  ///
  /// @param localCovariance Coverance matrix in local coordinates
  /// @param sigmaMax  Significance for the compatibility test
  BoundaryCheck(const ActsSymMatrixD<2>& localCovariance, double sigmaMax = 1);

  operator bool() const { return (m_type != Type::eNone); }
  bool operator!() const { return !bool(*this); }

  /// Check if the point is inside a polygon.
  ///
  /// @param point    Test point
  /// @param vertices Forward iterable container of convex polygon vertices.
  ///                 Calling `std::begin`/ `std::end` on the container must
  ///                 return an iterator where `*it` must be convertible to
  ///                 an `Acts::Vector2D`.
  ///
  /// The check takes into account whether tolerances or covariances are defined
  /// for the boundary check.
  template <typename Vector2DContainer>
  bool
  isInside(const Vector2D& point, const Vector2DContainer& vertices) const;

  /// Check if the point is inside a box aligned with the local axes.
  ///
  /// @param point   Test point
  /// @param loc0Min Lower bound along first axis
  /// @param loc0Max Upper bound along first axis
  /// @param loc1Min Lower bound along second axis
  /// @param loc1Max Upper bound along second axis
  ///
  /// The check takes into account whether tolerances or covariances are defined
  /// for the boundary check.
  bool
  isInside(const Vector2D& point,
           const Vector2D& lowerLeft,
           const Vector2D& upperRight) const;

  /// Calculate the signed, weighted, closest distance to a polygonal boundary.
  ///
  /// @param point Test point
  /// @param vertices Forward iterable container of convex polygon vertices.
  ///                 Calling `std::begin`/ `std::end` on the container must
  ///                 return an iterator where `*it` must be convertible to
  ///                 an `Acts::Vector2D`.
  /// @return Negative value if inside, positive if outside
  ///
  /// If a covariance is defined, the distance is the corresponding Mahalanobis
  /// distance. Otherwise, it is the Eucleadian distance.
  template <typename Vector2DContainer>
  double
  distance(const Vector2D& point, const Vector2DContainer& vertices) const;

  /// Calculate the signed, weighted, closest distance to an aligned box.
  ///
  /// @param point Test point
  /// @param loc0Min Minimal value along the first local axis
  /// @param loc0Max Maximal value along the first local axis
  /// @param loc1Min Minimal value along the first local axis
  /// @param loc1Max Maximal value along the first local axis
  /// @return Negative value if inside, positive if outside
  ///
  /// If a covariance is defined, the distance is the corresponding Mahalanobis
  /// distance. Otherwise, it is the Eucleadian distance.
  double
  distance(const Vector2D& point,
           const Vector2D& lowerLeft,
           const Vector2D& upperRight) const;

private:
  enum class Type {
    eNone,      ///< disable boundary check
    eAbsolute,  ///< absolute cut
    eChi2       ///< chi2-based cut with full correlations
  };

  /// Return a new BoundaryCheck with updated covariance.
  /// @param jacobian Tranform Jacobian for the covariance
  /// @warning This currently only transforms the covariance and does not work
  ///          for the tolerance based check.
  BoundaryCheck
  transformed(const ActsMatrixD<2, 2>& jacobian) const;

  /// Check if the point is inside the polygon w/o any tolerances.
  template <typename Vector2DContainer>
  bool
  isInsidePolygon(const Vector2D&          point,
                  const Vector2DContainer& vertices) const;

  /// Check if the point is inside the aligned box
  bool
  isInsideRectangle(const Vector2D& point,
                    const Vector2D& lowerLeft,
                    const Vector2D& upperRight) const;

  /// Check if the distance vector is within the absolute or relative limits.
  bool
  isTolerated(const Vector2D& delta) const;

  /// Compute vector norm based on the covariance.
  double
  squaredNorm(const Vector2D& x) const;

  /// Calculate the closest point on the polygon.
  template <typename Vector2DContainer>
  Vector2D
  computeClosestPointOnPolygon(const Vector2D&          point,
                               const Vector2DContainer& vertices) const;

  /// Calculate the closest point on the box
  Vector2D
  computeEuclideanClosestPointOnRectangle(const Vector2D& point,
                                          const Vector2D& lowerLeft,
                                          const Vector2D& upperRight) const;

  /// metric weight matrix: identity for absolute mode or inverse covariance
  ActsSymMatrixD<2> m_weight;

  /// dual use: absolute tolerances or relative chi2/ sigma cut.
  Vector2D m_tolerance;
  Type     m_type;

  // To acces the m_type
  friend class CylinderBounds;
  friend class RectangleBounds;
  // To be able to use `transformed`
  friend class CylinderBounds;
  friend class DiscTrapezoidalBounds;
  // EllipseBounds needs a custom implementation
  friend class EllipseBounds;
};

}  // namespace Acts

inline Acts::BoundaryCheck::BoundaryCheck(bool check)
  : m_weight(ActsSymMatrixD<2>::Identity())
  , m_tolerance(0, 0)
  , m_type(check ? Type::eAbsolute : Type::eNone)
{
}

inline Acts::BoundaryCheck::BoundaryCheck(bool   checkLocal0,
                                          bool   checkLocal1,
                                          double tolerance0,
                                          double tolerance1)
  : m_weight(ActsSymMatrixD<2>::Identity())
  , m_tolerance(checkLocal0 ? tolerance0 : DBL_MAX,
                checkLocal1 ? tolerance1 : DBL_MAX)
  , m_type(Type::eAbsolute)
{
}

inline Acts::BoundaryCheck::BoundaryCheck(
    const ActsSymMatrixD<2>& localCovariance,
    double                   sigmaMax)
  : m_weight(localCovariance.inverse())
  , m_tolerance(sigmaMax, 0)
  , m_type(Type::eChi2)
{
}

inline Acts::BoundaryCheck
Acts::BoundaryCheck::transformed(const ActsMatrixD<2, 2>& jacobian) const
{
  BoundaryCheck bc = *this;
  if (m_type == Type::eAbsolute) {
    // project tolerances to the new system. depending on the jacobian we need
    // to check both tolerances, even when the initial check does not.
    bc.m_tolerance = (jacobian * m_tolerance).cwiseAbs();
  } else /* Type::eChi2 */ {
    bc.m_weight
        = (jacobian * m_weight.inverse() * jacobian.transpose()).inverse();
  }
  return bc;
}

template <typename Vector2DContainer>
inline bool
Acts::BoundaryCheck::isInside(const Vector2D&          point,
                              const Vector2DContainer& vertices) const
{
  // a compatible point must be either completely in the polygon or on the
  // outside but within the defined tolerance relative to the closest point
  if (isInsidePolygon(point, vertices)) {
    return true;
  } else {
    auto closestPoint = computeClosestPointOnPolygon(point, vertices);
    return isTolerated(closestPoint - point);
  }
}

template <typename Vector2DContainer>
inline bool
Acts::BoundaryCheck::isInsidePolygon(const Vector2D&          point,
                                     const Vector2DContainer& vertices) const
{
  // when we move along the edges of a convex polygon, a point on the inside of
  // the polygon will always appear on the same side of each edge.
  // a point on the outside will switch sides at least once.

  // returns which side of the connecting line between `ll0` and `ll1` the point
  // `p` is on. computes the sign of the z-component of the cross-product
  // between the line normal vector and the vector from `ll0` to `p`.
  auto lineSide = [&](auto&& ll0, auto&& ll1) {
    auto normal = ll1 - ll0;
    auto delta  = point - ll0;
    return std::signbit((normal[0] * delta[1]) - (normal[1] * delta[0]));
  };

  auto     iv = std::begin(vertices);
  Vector2D l0 = *iv;
  Vector2D l1 = *(++iv);
  // use vertex0 to vertex1 to define reference sign and compare w/ all edges
  auto reference = lineSide(l0, l1);
  for (++iv; iv != std::end(vertices); ++iv) {
    l0 = l1;
    l1 = *iv;
    if (lineSide(l0, l1) != reference) {
      return false;
    }
  }
  // manual check for last edge from last vertex back to the first vertex
  if (lineSide(l1, *std::begin(vertices)) != reference) {
    return false;
  }
  // point was always on the same side. point must be inside.
  return true;
}

inline bool
Acts::BoundaryCheck::isInsideRectangle(const Vector2D& point,
                                       const Vector2D& lowerLeft,
                                       const Vector2D& upperRight) const
{
  return (lowerLeft[0] <= point[0]) && (point[0] < upperRight[0])
      && (lowerLeft[1] <= point[1]) && (point[1] < upperRight[1]);
}

inline bool
Acts::BoundaryCheck::isInside(const Vector2D& point,
                              const Vector2D& lowerLeft,
                              const Vector2D& upperRight) const
{
  if (isInsideRectangle(point, lowerLeft, upperRight)) {
    return true;
  } else {

    Vector2D closestPoint;

    if (m_type == Type::eNone || m_type == Type::eAbsolute) {
      // absolute, can calculate directly
      closestPoint = computeEuclideanClosestPointOnRectangle(
          point, lowerLeft, upperRight);

    } else /* Type::eChi2 */ {
      // need to calculate by projection and squarednorm
      Vector2D vertices[] = {{lowerLeft[0], lowerLeft[1]},
                             {upperRight[0], lowerLeft[1]},
                             {upperRight[0], upperRight[1]},
                             {lowerLeft[0], upperRight[1]}};
      closestPoint = computeClosestPointOnPolygon(point, vertices);
    }

    return isTolerated(closestPoint - point);
  }
}

template <typename Vector2DContainer>
inline double
Acts::BoundaryCheck::distance(const Acts::Vector2D&    point,
                              const Vector2DContainer& vertices) const
{
  // TODO 2017-04-06 msmk: this should be calculable directly
  double d = squaredNorm(point - computeClosestPointOnPolygon(point, vertices));
  d        = std::sqrt(d);
  return isInsidePolygon(point, vertices) ? -d : d;
}

inline double
Acts::BoundaryCheck::distance(const Acts::Vector2D& point,
                              const Vector2D&       lowerLeft,
                              const Vector2D&       upperRight) const
{
  if (m_type == Type::eNone || m_type == Type::eAbsolute) {

    // compute closest point on box
    double d = (point - computeEuclideanClosestPointOnRectangle(
                            point, lowerLeft, upperRight))
                   .norm();
    return isInsideRectangle(point, lowerLeft, upperRight) ? -d : d;

  } else /* Type::eChi2 */ {
    Vector2D vertices[] = {{lowerLeft[0], lowerLeft[1]},
                           {upperRight[0], lowerLeft[1]},
                           {upperRight[0], upperRight[1]},
                           {lowerLeft[0], upperRight[1]}};
    return distance(point, vertices);
  }
}

inline bool
Acts::BoundaryCheck::isTolerated(const Vector2D& delta) const
{
  if (m_type == Type::eNone) {
    return true;
  } else if (m_type == Type::eAbsolute) {
    return (std::abs(delta[0]) <= m_tolerance[0])
        && (std::abs(delta[1]) <= m_tolerance[1]);
  } else /* Type::eChi2 */ {
    // Mahalanobis distances mean is 2 in 2-dim. cut is 1-d sigma.
    return (squaredNorm(delta) < (2 * m_tolerance[0]));
  }
}

inline double
Acts::BoundaryCheck::squaredNorm(const Vector2D& x) const
{
  return (x.transpose() * m_weight * x).value();
}

template <typename Vector2DContainer>
inline Acts::Vector2D
Acts::BoundaryCheck::computeClosestPointOnPolygon(
    const Acts::Vector2D&    point,
    const Vector2DContainer& vertices) const
{
  // calculate the closest position on the segment between `ll0` and `ll1` to
  // the point as measured by the metric induced by the weight matrix
  auto closestOnSegment = [&](auto&& ll0, auto&& ll1) {
    // normal vector and position of the closest point along the normal
    auto n = ll1 - ll0;
    auto f = (n.transpose() * m_weight * n).value();
    auto u = std::isnormal(f)
        ? -((ll0 - point).transpose() * m_weight * n).value() / f
        : 0.5;  // ll0 and ll1 are so close it doesn't matter
    // u must be in [0, 1] to still be on the polygon segment
    return ll0 + std::min(std::max(u, 0.0), 1.0) * n;
  };

  auto     iv      = std::begin(vertices);
  Vector2D l0      = *iv;
  Vector2D l1      = *(++iv);
  Vector2D closest = closestOnSegment(l0, l1);
  // Calculate the closest point on other connecting lines and compare distances
  for (++iv; iv != std::end(vertices); ++iv) {
    l0               = l1;
    l1               = *iv;
    Vector2D current = closestOnSegment(l0, l1);
    if (squaredNorm(current - point) < squaredNorm(closest - point)) {
      closest = current;
    }
  }
  // final edge from last vertex back to the first vertex
  Vector2D last = closestOnSegment(l1, *std::begin(vertices));
  if (squaredNorm(last - point) < squaredNorm(closest - point)) {
    closest = last;
  }
  return closest;
}

inline Acts::Vector2D
Acts::BoundaryCheck::computeEuclideanClosestPointOnRectangle(
    const Vector2D& point,
    const Vector2D& lowerLeft,
    const Vector2D& upperRight) const
{

  /*
   *
   *        |                 |
   *   IV   |       V         | I
   *        |                 |
   *  ------------------------------
   *        |                 |
   *        |                 |
   *   VIII |     INSIDE      | VI
   *        |                 |
   *        |                 |
   *  ------------------------------
   *        |                 |
   *   III  |      VII        | II
   *        |                 |
   *
   */

  double l0 = point[0], l1 = point[1];
  double loc0Min = lowerLeft[0], loc0Max = upperRight[0];
  double loc1Min = lowerLeft[1], loc1Max = upperRight[1];

  // check if inside
  if (loc0Min <= l0 && l0 < loc0Max && loc1Min <= l1 && l1 < loc1Max) {
    // INSIDE
    double   dist = std::abs(loc0Max - l0);
    Vector2D cls(loc0Max, l1);

    double test = std::abs(loc0Min - l0);
    if (test <= dist) {
      dist = test;
      cls  = {loc0Min, l1};
    }

    test = std::abs(loc1Max - l1);
    if (test <= dist) {
      dist = test;
      cls  = {l0, loc1Max};
    }

    test = std::abs(loc1Min - l1);
    if (test <= dist) {
      return {l0, loc1Min};
    }
    return cls;
  } else {
    // OUTSIDE, check sectors
    if (l0 > loc0Max) {
      if (l1 > loc1Max) {  // I
        return {loc0Max, loc1Max};
      } else if (l1 <= loc1Min) {  // II
        return {loc0Max, loc1Min};
      } else {  // VI
        return {loc0Max, l1};
      }
    } else if (l0 < loc0Min) {
      if (l1 > loc1Max) {  // IV
        return {loc0Min, loc1Max};
      } else if (l1 <= loc1Min) {  // III
        return {loc0Min, loc1Min};
      } else {  // VIII
        return {loc0Min, l1};
      }
    } else {
      if (l1 > loc1Max) {  // V
        return {l0, loc1Max};
      } else {  // l1 <= loc1Min # VII
        return {l0, loc1Min};
      }
      // third case not necessary, see INSIDE above
    }
  }
}
