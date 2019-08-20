// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorElementBase.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Plugins/Identification/Identifier.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
class Surface;
class PlanarBounds;
class DiscBounds;
class SurfaceMaterial;
class DigitizationModule;
}

namespace FW {

namespace Generic {

  /// @class GenericDetectorElement
  ///
  /// This is a lightweight type of detector element,
  /// it simply implements the base class.
  ///
  class GenericDetectorElement : public Acts::IdentifiedDetectorElement
  {
  public:
    /// Constructor for single sided detector element
    /// - bound to a Plane Surface
    ///
    /// @param identifier is the module identifier
    /// @param transform is the transform that element the layer in 3D frame
    /// @param pBounds is the planar bounds for the planar detector element
    /// @param thickness is the module thickness
    /// @param material is the (optional) Surface material associated to it
    GenericDetectorElement(
        const Identifier                                identifier,
        std::shared_ptr<const Acts::Transform3D>        transform,
        std::shared_ptr<const Acts::PlanarBounds>       pBounds,
        double                                          thickness,
        std::shared_ptr<const Acts::SurfaceMaterial>    material = nullptr,
        std::shared_ptr<const Acts::DigitizationModule> digitzationModule
        = nullptr);

    /// Constructor for single sided detector element
    /// - bound to a Disc Surface
    ///
    /// @param identifier is the module identifier
    /// @param transform is the transform that element the layer in 3D frame
    /// @param dBounds is the planar bounds for the disc like detector element
    /// @param thickness is the module thickness
    /// @param material is the (optional) Surface material associated to it
    GenericDetectorElement(
        const Identifier                                identifier,
        std::shared_ptr<const Acts::Transform3D>        transform,
        std::shared_ptr<const Acts::DiscBounds>         dBounds,
        double                                          thickness,
        std::shared_ptr<const Acts::SurfaceMaterial>    material = nullptr,
        std::shared_ptr<const Acts::DigitizationModule> digitzationModule
        = nullptr);

    /// Identifier
    Identifier
    identifier() const override final;

    /// Return local to global transform associated with this identifier
    ///
    /// @note this is called from the surface().transform() in the PROXY mode
    const Acts::Transform3D&
    transform() const final override;

    /// Return surface associated with this identifier,
    const Acts::Surface&
    surface() const final override;

    /// Set the identifier after construction (sometimes needed)
    void
    assignIdentifier(const Identifier& identifier);

    /// The maximal thickness of the detector element wrt normal axis
    double
    thickness() const final override;

    /// Retrieve the DigitizationModule
    const std::shared_ptr<const Acts::DigitizationModule>
    digitizationModule() const final override;

  private:
    /// the element representation
    /// identifier
    Identifier m_elementIdentifier;
    /// the transform for positioning in 3D space
    std::shared_ptr<const Acts::Transform3D> m_elementTransform;
    /// the surface represented by it
    std::shared_ptr<const Acts::Surface> m_elementSurface;
    /// the element thickness
    double m_elementThickness;
    /// store either
    std::shared_ptr<const Acts::PlanarBounds> m_elementPlanarBounds = nullptr;
    std::shared_ptr<const Acts::DiscBounds>   m_elementDiscBounds   = nullptr;
    /// The Digitization module
    std::shared_ptr<const Acts::DigitizationModule> m_digitizationModule
        = nullptr;
  };

  inline void
  FW::Generic::GenericDetectorElement::assignIdentifier(
      const Identifier& identifier)
  {
    m_elementIdentifier = identifier;
  }

  inline Identifier
  FW::Generic::GenericDetectorElement::identifier() const
  {
    return m_elementIdentifier;
  }

  inline const Acts::Transform3D&
  FW::Generic::GenericDetectorElement::transform() const
  {
    return *m_elementTransform;
  }

  inline const Acts::Surface&
  FW::Generic::GenericDetectorElement::surface() const
  {
    return *m_elementSurface;
  }

  inline double
  FW::Generic::GenericDetectorElement::thickness() const
  {
    return m_elementThickness;
  }

  inline const std::shared_ptr<const Acts::DigitizationModule>
  FW::Generic::GenericDetectorElement::digitizationModule() const
  {
    return m_digitizationModule;
  }

}  // end of namespace Generic

}  // end of namespace FW
