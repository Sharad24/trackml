// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <vector>
#include "ACTFW/Barcode/Barcode.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

namespace FW {

/// Typedef the pdg code
typedef int pdg_type;

namespace Data {

  /// @brief Particle information struct for physics process samplers:
  /// - all quatities are calculated at first construction as they may
  ///   be used by downstream samplers
  ///
  /// @note if a sampler changes one of the parameters, consistency
  /// can be broken, so it should update the rest (no checking done)
  class SimParticle
  {

  public:
    /// @brief Default Constructor
    SimParticle() = default;

    /// @brief Construct a particle consistently
    ///
    /// @param position The particle position at construction
    /// @param momentum The particle momentum at construction
    /// @param m The particle mass
    /// @param q The partilce charge
    /// @param barcode The particle barcode
    /// @param tStamp is the current time stamp
    SimParticle(const Acts::Vector3D& position,
                const Acts::Vector3D& momentum,
                double                m,
                double                q,
                pdg_type              pdg     = 0,
                barcode_type          barcode = 0,
                double                tStamp  = 0.)
      : m_position(position)
      , m_momentum(momentum)
      , m_m(m)
      , m_q(q)
      , m_p(momentum.norm())
      , m_pT(Acts::VectorHelpers::perp(momentum))
      , m_pdg(pdg)
      , m_barcode(barcode)
      , m_timeStamp(tStamp)
    {
      m_E     = std::sqrt(m_p * m_p + m_m * m_m);
      m_beta  = (m_p / m_E);
      m_gamma = (m_E / m_m);
    }

    /// Default
    SimParticle(const SimParticle& sp) = default;

    /// @brief Set the limits
    ///
    /// @param x0Limit the limit in X0 to be passed
    /// @param l0Limit the limit in L0 to be passed
    /// @param timeLimit the readout time limit to be passed
    void
    setLimits(double x0Limit,
              double l0Limit,
              double timeLimit = std::numeric_limits<double>::max())
    {
      m_limitInX0 = x0Limit;
      m_limitInL0 = l0Limit;
      m_timeLimit = timeLimit;
    }

    /// @brief Place the particle int he detector and set barcode
    ///
    /// @param deltaE is the energy loss to be applied
    void
    place(Acts::Vector3D position, barcode_type barcode, double timeStamp = 0.)
    {
      m_position  = std::move(position);
      m_barcode   = barcode;
      m_timeStamp = timeStamp;
    }

    /// @brief Update the particle with applying energy loss
    ///
    /// @param deltaE is the energy loss to be applied
    void
    scatter(Acts::Vector3D nmomentum)
    {
      m_momentum = std::move(nmomentum);
      m_pT       = Acts::VectorHelpers::perp(m_momentum);
    }

    /// @brief Update the particle with applying energy loss
    ///
    /// @param deltaE is the energy loss to be applied
    void
    energyLoss(double deltaE)
    {
      // particle falls to rest
      if (m_E - deltaE < m_m) {
        m_E        = m_m;
        m_p        = 0.;
        m_pT       = 0.;
        m_beta     = 0.;
        m_gamma    = 1.;
        m_momentum = Acts::Vector3D(0., 0., 0.);
        m_alive    = false;
      }
      // updatet the parameters
      m_E -= deltaE;
      m_p        = std::sqrt(m_E * m_E - m_m * m_m);
      m_momentum = m_p * m_momentum.normalized();
      m_pT       = Acts::VectorHelpers::perp(m_momentum);
      m_beta     = (m_p / m_E);
      m_gamma    = (m_E / m_m);
    }

    /// @brief Update the particle with a new position and momentum,
    /// this corresponds to a step update
    ///
    /// @param position New position after update
    /// @param momentum New momentum after update
    /// @param deltaPathX0 passed since last step
    /// @param deltaPathL0 passed since last step
    /// @param deltaTime The time elapsed
    ///
    /// @return break condition
    bool
    update(const Acts::Vector3D& position,
           const Acts::Vector3D& momentum,
           double                deltaPathX0 = 0.,
           double                deltaPathL0 = 0.,
           double                deltaTime   = 0.)
    {
      m_position = position;
      m_momentum = momentum;
      m_p        = momentum.norm();
      if (m_p) {
        m_pT = Acts::VectorHelpers::perp(momentum);
        m_E  = std::sqrt(m_p * m_p + m_m * m_m);
        m_timeStamp += deltaTime;
        m_beta  = (m_p / m_E);
        m_gamma = (m_E / m_m);

        // set parameters and check limits
        m_pathInX0 += deltaPathX0;
        m_pathInL0 += deltaPathL0;
        m_timeStamp += deltaTime;
        if (m_pathInX0 >= m_limitInX0 || m_pathInL0 >= m_limitInL0
            || m_timeStamp > m_timeLimit) {
          m_alive = false;
        }
      }
      return !m_alive;
    }

    /// @bref boost the particle
    // void boost(){
    //
    // }

    /// @brief Access methods: position
    const Acts::Vector3D&
    position() const
    {
      return m_position;
    }

    /// @brief Access methods: momentum
    const Acts::Vector3D&
    momentum() const
    {
      return m_momentum;
    }

    /// @brief Access methods: p
    const double
    p() const
    {
      return m_p;
    }

    /// @brief Access methods: pT
    const double
    pT() const
    {
      return m_pT;
    }

    /// @brief Access methods: E
    const double
    E() const
    {
      return m_E;
    }

    /// @brief Access methods: m
    const double
    m() const
    {
      return m_m;
    }

    /// @brief Access methods: beta
    const double
    beta() const
    {
      return m_beta;
    }

    /// @brief Access methods: gamma
    const double
    gamma() const
    {
      return m_gamma;
    }

    /// @brief Access methods: charge
    const double
    q() const
    {
      return m_q;
    }

    /// @brief Access methods: pdg code
    const pdg_type
    pdg() const
    {
      return m_pdg;
    }

    /// @brief Access methods: barcode
    const barcode_type
    barcode() const
    {
      return m_barcode;
    }

    /// @brief Access methods: path/X0
    const double
    pathInX0() const
    {
      return m_pathInX0;
    }

    /// @brief Access methods: limit/X0
    const double
    limitInX0() const
    {
      return m_limitInX0;
    }

    /// @brief Access methods: pdg code
    const double
    pathInL0() const
    {
      return m_limitInX0;
    }

    /// @brief Access methods: barcode
    const double
    limitInL0() const
    {
      return m_limitInL0;
    }

    /// @brief boolean operator indicating the particle to be alive
    operator bool() { return m_alive; }

  private:
    Acts::Vector3D m_position = Acts::Vector3D(0., 0., 0.);  //!< kinematic info
    Acts::Vector3D m_momentum = Acts::Vector3D(0., 0., 0.);  //!< kinematic info

    double       m_m       = 0.;  //!< particle mass
    double       m_E       = 0.;  //!< total energy
    double       m_q       = 0.;  //!< the charge
    double       m_beta    = 0.;  //!< relativistic beta factor
    double       m_gamma   = 1.;  //!< relativistic gamma factor
    double       m_p       = 0.;  //!< momentum magnitude
    double       m_pT      = 0.;  //!< transverse momentum magnitude
    pdg_type     m_pdg     = 0;   //!< pdg code of the particle
    barcode_type m_barcode = 0;   //!< barcode of the particle

    double m_pathInX0 = 0.;  //!< passed path in X0
    double m_limitInX0
        = std::numeric_limits<double>::max();  //!< path limit in X0

    double m_pathInL0 = 0.;  //!< passed path in L0
    double m_limitInL0
        = std::numeric_limits<double>::max();  //!< path limit in X0

    double m_timeStamp = 0.;  //!< passed time elapsed
    double m_timeLimit = std::numeric_limits<double>::max();  // time limit

    bool m_alive = true;  //!< the particle is alive
  };

}  // end of namespace Data
}  // end of namespace FW
