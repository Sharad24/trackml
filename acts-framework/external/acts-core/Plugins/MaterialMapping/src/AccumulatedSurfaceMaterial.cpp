// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AccumulatedSurfaceMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialMapping/AccumulatedSurfaceMaterial.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"

// Default Constructor - for homogeneous material
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(double splitFactor)
  : m_splitFactor(splitFactor)
{
  AccumulatedVector accMat = {{AccumulatedMaterialProperties()}};
  m_accumulatedMaterial    = {{accMat}};
}

// Binned Material constructor with split factor
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(
    const BinUtility& binUtility,
    double            splitFactor)
  : m_binUtility(binUtility), m_splitFactor(splitFactor)
{
  size_t            bins0 = m_binUtility.bins(0);
  size_t            bins1 = m_binUtility.bins(1);
  AccumulatedVector accVec(bins0, AccumulatedMaterialProperties());
  m_accumulatedMaterial = AccumulatedMatrix(bins1, accVec);
}

// Assign a material properites object
void
Acts::AccumulatedSurfaceMaterial::accumulate(const Vector2D&           lp,
                                             const MaterialProperties& mp,
                                             double pathCorrection)
{
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
  } else {
    size_t bin0 = m_binUtility.bin(lp, 0);
    size_t bin1 = m_binUtility.bin(lp, 1);
    m_accumulatedMaterial[bin1][bin0].accumulate(mp, pathCorrection);
  }
}

// Assign a material properites object
void
Acts::AccumulatedSurfaceMaterial::accumulate(const Vector3D&           gp,
                                             const MaterialProperties& mp,
                                             double pathCorrection)
{
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
  } else {
    std::array<size_t, 3> bTriple = m_binUtility.binTriple(gp);
    m_accumulatedMaterial[bTriple[1]][bTriple[0]].accumulate(mp,
                                                             pathCorrection);
  }
}

// Assign a vector of material properites object
void
Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector3D& gp,
    const std::vector<std::pair<MaterialProperties, Vector3D>>& mps,
    double pathCorrection)
{
  if (m_binUtility.dimensions() == 0) {
    for (auto mp : mps) {
      m_accumulatedMaterial[0][0].accumulate(mp.first, pathCorrection);
    }
  } else {
    std::array<size_t, 3> bTriple = m_binUtility.binTriple(gp);
    for (auto mp : mps) {
      m_accumulatedMaterial[bTriple[1]][bTriple[0]].accumulate(mp.first,
                                                               pathCorrection);
    }
  }
}

// Average the information accumulated during one event
void
Acts::AccumulatedSurfaceMaterial::eventAverage()
{
  for (auto& matVec : m_accumulatedMaterial) {
    for (auto& mat : matVec) {
      mat.eventAverage();
    }
  }
}

/// Total average creates SurfaceMaterial
std::unique_ptr<const Acts::SurfaceMaterial>
Acts::AccumulatedSurfaceMaterial::totalAverage()
{
  if (m_binUtility.bins() == 1) {
    // Return HomogeneousSurfaceMaterial
    return std::make_unique<HomogeneousSurfaceMaterial>(
        m_accumulatedMaterial[0][0].totalAverage().first, m_splitFactor);
  }
  // Create the properties matrix
  MaterialPropertiesMatrix mpMatrix(
      m_binUtility.bins(1),
      MaterialPropertiesVector(m_binUtility.bins(0), MaterialProperties()));
  // Loop over and fill
  for (size_t ib1 = 0; ib1 < m_binUtility.bins(1); ++ib1) {
    for (size_t ib0 = 0; ib0 < m_binUtility.bins(0); ++ib0) {
      mpMatrix[ib1][ib0] = m_accumulatedMaterial[ib1][ib0].totalAverage().first;
    }
  }
  // Now return the BinnedSurfaceMaterial
  return std::make_unique<const BinnedSurfaceMaterial>(
      m_binUtility, std::move(mpMatrix), m_splitFactor);
}
