// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace FW {

namespace Data {

  /// These are the SimParticles
  class SimParticle;

  /// @class SimIdentifier
  ///
  /// Identifier implementation for the ACTS framework
  /// including MC truth information
  class SimIdentifier
  {
  public:
    typedef unsigned long long identifier_type;
    typedef long long          identifier_diff;

    /// Default constructor
    SimIdentifier() = default;

    /// Constructor from identifier_type
    ///
    /// @param value is the identifier value
    explicit SimIdentifier(identifier_type value);

    /// Constructor from identifier_type
    ///
    /// @param value is the identifier value
    SimIdentifier(identifier_type                 value,
                  std::vector<const SimParticle*> truthParticles);

    /// Copy constructor
    ///
    /// @param other is the source identifier
    SimIdentifier(const SimIdentifier& other) = default;

    /// @param old is the assigment parameter
    SimIdentifier&
    operator=(const SimIdentifier& old)
        = default;

    /// @param value is the assigment parameter
    SimIdentifier&
    operator=(identifier_type value);

    /// Cast operators to value @todo to bool
    operator identifier_type() const { return m_id; }
    identifier_type
    value() const
    {
      return m_id;
    }

    /// Attach a truth particle
    ///
    /// @param particle is the truth particle to be attached;
    void
    attachTruthParticle(const SimParticle* particle);

    /// Read the truth particles
    const std::vector<const SimParticle*>&
    truthParticles() const;

    /// @param other is the comparison parameter
    bool
    operator==(const SimIdentifier& other) const;

    /// @param other is the comparison parameter
    bool
    operator!=(const SimIdentifier& other) const;

    /// @param other is the comparison parameter
    bool
    operator<(const SimIdentifier& other) const;

    /// @param other is the comparison parameter
    bool
    operator>(const SimIdentifier& other) const;

    /// @param other is the comparison parameter
    bool
    operator<=(const SimIdentifier& other) const;

    /// @param other is the comparison parameter
    bool
    operator>=(const SimIdentifier& other) const;

  private:
    identifier_type                 m_id = 0;  //! the store identifier value
    std::vector<const SimParticle*> m_truthParticles
        = {};  //!< the attached particles
  };

  inline SimIdentifier::SimIdentifier(identifier_type value)
    : m_id(value), m_truthParticles()
  {
  }

  inline SimIdentifier::SimIdentifier(
      identifier_type                 value,
      std::vector<const SimParticle*> truthParticles)
    : m_id(value), m_truthParticles(std::move(truthParticles))
  {
  }

  inline SimIdentifier&
  SimIdentifier::operator=(identifier_type value)
  {
    m_id = value;
    return (*this);
  }

  inline bool
  SimIdentifier::operator==(const SimIdentifier& other) const
  {
    return (m_id == other.m_id);
  }

  inline bool
  SimIdentifier::operator!=(const SimIdentifier& other) const
  {
    return (m_id != other.m_id);
  }

  inline bool
  SimIdentifier::operator<(const SimIdentifier& other) const
  {
    return (m_id < other.m_id);
  }

  inline bool
  SimIdentifier::operator>(const SimIdentifier& other) const
  {
    return (m_id > other.m_id);
  }

  inline bool
  SimIdentifier::operator<=(const SimIdentifier& other) const
  {
    return (m_id <= other.m_id);
  }

  inline bool
  SimIdentifier::operator>=(const SimIdentifier& other) const
  {
    return (m_id >= other.m_id);
  }

  inline void
  SimIdentifier::attachTruthParticle(const SimParticle* particle)
  {
    m_truthParticles.push_back(particle);
  }

  inline const std::vector<const SimParticle*>&
  SimIdentifier::truthParticles() const
  {
    return m_truthParticles;
  }

}  // end of namespace Test
}  // end of namespace FW

using Identifier = FW::Data::SimIdentifier;
