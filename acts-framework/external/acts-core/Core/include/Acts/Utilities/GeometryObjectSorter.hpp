// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeometryObjectSorter.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Definitions.hpp"

namespace Acts {

// @class ObjectSorterT
///
template <class T>
class ObjectSorterT : public std::binary_function<T, T, bool>
{
public:
  /// Constructor from a binning value
  ///
  /// @param bValue is the value in which the binning is done
  /// @param transform is an optional transform to be performed
  ObjectSorterT(BinningValue bValue) : m_binningValue(bValue) {}

  /// Comparison operator
  ///
  /// @tparam one first object
  /// @tparam two second object
  ///
  /// @return boolen indicator
  bool
  operator()(T one, T two) const
  {
    using Acts::VectorHelpers::perp;
    using Acts::VectorHelpers::phi;
    using Acts::VectorHelpers::eta;
    // switch the binning value
    // - binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
    switch (m_binningValue) {
    // compare on x
    case binX: {
      return (one.x() < two.x());
    }
    // compare on y
    case binY: {
      return (one.y() < two.y());
    }
    // compare on z
    case binZ: {
      return (one.z() < two.z());
    }
    // compare on r
    case binR: {
      return (perp(one) < perp(two));
    }
    // compare on phi
    case binPhi: {
      return (phi(one) < phi(two));
    }
    // compare on eta
    case binEta: {
      return (eta(one) < eta(two));
    }
    // default for the moment
    default: {
      return (one.norm() < two.norm());
    }
    }
  }

  BinningValue
  binningValue() const
  {
    return m_binningValue;
  }

private:
  BinningValue m_binningValue;  ///< the binning value
};

/// @class DistanceSorterT
///
/// This will check on absolute distance
template <class T>
class DistanceSorterT : public std::binary_function<T, T, bool>
{
public:
  /// Constructor from a binning value
  ///
  /// @param bValue is the value in which the binning is done
  /// @param reference is the reference point
  DistanceSorterT(BinningValue bValue, Vector3D reference)
    : m_binningValue(bValue)
    , m_reference(reference)
    , m_refR(VectorHelpers::perp(reference))
    , m_refPhi(VectorHelpers::phi(reference))
    , m_refEta(VectorHelpers::eta(reference))
  {
  }

  /// Comparison operator
  ///
  /// @tparam one first object
  /// @tparam two second object
  ///
  /// @return boolen indicator
  bool
  operator()(T one, T two) const
  {
    using Acts::VectorHelpers::perp;
    using Acts::VectorHelpers::phi;
    using Acts::VectorHelpers::eta;
    // switch the binning value
    // - binX, binY, binZ, binR, binPhi, binRPhi, binH, binEta
    switch (m_binningValue) {
    // compare on diff x
    case binX: {
      double diffOneX = one.x() - m_reference.x();
      double diffTwoX = two.x() - m_reference.x();
      return (diffOneX * diffOneX < diffTwoX * diffTwoX);
    }
    // compare on diff y
    case binY: {
      double diffOneY = one.y() - m_reference.y();
      double diffTwoY = two.y() - m_reference.y();
      return (diffOneY * diffOneY < diffTwoY * diffTwoY);
    }
    // compare on diff z
    case binZ: {
      double diffOneZ = one.z() - m_reference.z();
      double diffTwoZ = two.z() - m_reference.z();
      return (diffOneZ * diffOneZ < diffTwoZ * diffTwoZ);
    }
    // compare on r
    case binR: {
      double diffOneR = perp(one) - m_refR;
      double diffTwoR = perp(two) - m_refR;
      return (diffOneR * diffOneR < diffTwoR * diffTwoR);
    }
    // compare on phi /// @todo add cyclic value
    case binPhi: {
      double diffOnePhi = phi(one) - m_refPhi;
      double diffTwoPhi = phi(two) - m_refPhi;
      return (diffOnePhi * diffOnePhi < diffTwoPhi * diffTwoPhi);
    }
    // compare on eta
    case binEta: {
      double diffOneEta = eta(one) - m_refEta;
      double diffTwoEta = eta(two) - m_refEta;
      return (diffOneEta * diffOneEta < diffTwoEta * diffTwoEta);
    }
    // default for the moment
    default: {
      T diffOne(one - m_reference);
      T diffTwo(two - m_reference);
      return (diffOne.mag2() < diffTwo.mag2());
    }
    }
  }

private:
  BinningValue m_binningValue;  ///< the binning value
  T            m_reference;
  double       m_refR;
  double       m_refPhi;
  double       m_refEta;
};

/// @class GeometryObjectSorter
///
template <class T>
class GeometryObjectSorterT : public std::binary_function<T, T, bool>
{
public:
  /// Constructor from a binning value
  ///
  /// @param bValue is the value in which the binning is done
  /// @param transform is an optional transform to be performed
  GeometryObjectSorterT(BinningValue                       bValue,
                        std::shared_ptr<const Transform3D> transform = nullptr)
    : m_objectSorter(bValue), m_transform(std::move(transform))
  {
  }

  /// Comparison operator
  ///
  /// @tparam one first object
  /// @tparam two second object
  ///
  /// @return boolen indicator
  bool
  operator()(T one, T two) const
  {
    // get the pos one / pos two
    Vector3D posOne = m_transform
        ? m_transform->inverse()
            * one->binningPosition(m_objectSorter.binningValue())
        : one->binningPosition(m_objectSorter.binningValue());
    Vector3D posTwo = m_transform
        ? m_transform->inverse()
            * two->binningPosition(m_objectSorter.binningValue())
        : two->binningPosition(m_objectSorter.binningValue());
    // now call the distance sorter
    return m_objectSorter.operator()(posOne, posTwo);
  }

protected:
  ObjectSorterT<Vector3D>            m_objectSorter;
  std::shared_ptr<const Transform3D> m_transform;
};
}