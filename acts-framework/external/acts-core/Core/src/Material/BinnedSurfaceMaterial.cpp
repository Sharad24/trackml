// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinnedSurfaceMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const BinUtility&        binUtility,
    MaterialPropertiesVector fullProperties,
    double                   splitFactor)
  : SurfaceMaterial(splitFactor), m_binUtility(binUtility)
{
  // fill the material with deep copy
  m_fullMaterial.push_back(std::move(fullProperties));
}

Acts::BinnedSurfaceMaterial::BinnedSurfaceMaterial(
    const BinUtility&        binUtility,
    MaterialPropertiesMatrix fullProperties,
    double                   splitFactor)
  : Acts::SurfaceMaterial(splitFactor)
  , m_binUtility(binUtility)
  , m_fullMaterial(std::move(fullProperties))
{
}

Acts::BinnedSurfaceMaterial&
Acts::BinnedSurfaceMaterial::operator*=(double scale)
{
  for (auto& materialVector : m_fullMaterial) {
    for (auto& materialBin : materialVector) {
      (materialBin) *= scale;
    }
  }
  return (*this);
}

const Acts::MaterialProperties&
Acts::BinnedSurfaceMaterial::materialProperties(const Vector2D& lp) const
{
  // the first bin
  size_t ibin0 = m_binUtility.bin(lp, 0);
  size_t ibin1 = m_binUtility.max(1) != 0u ? m_binUtility.bin(lp, 1) : 0;
  return m_fullMaterial[ibin1][ibin0];
}

const Acts::MaterialProperties&
Acts::BinnedSurfaceMaterial::materialProperties(const Acts::Vector3D& gp) const
{
  // the first bin
  size_t ibin0 = m_binUtility.bin(gp, 0);
  size_t ibin1 = m_binUtility.max(1) != 0u ? m_binUtility.bin(gp, 1) : 0;
  return m_fullMaterial[ibin1][ibin0];
}

std::ostream&
Acts::BinnedSurfaceMaterial::dump(std::ostream& sl) const
{
  sl << "Acts::BinnedSurfaceMaterial : " << std::endl;
  sl << "   - Number of Material bins [0,1] : " << m_binUtility.max(0) + 1
     << " / " << m_binUtility.max(1) + 1 << std::endl;
  sl << "   - Parse full update material    : " << std::endl;  //
  // output  the full material
  unsigned int imat1 = 0;
  for (auto& materialVector : m_fullMaterial) {
    unsigned int imat0 = 0;
    // the vector iterator
    for (auto& materialBin : materialVector) {
      sl << " Bin [" << imat1 << "][" << imat0 << "] - " << (materialBin);
      ++imat0;
    }
    ++imat1;
  }
  return sl;
}
