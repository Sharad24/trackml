// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2016-05-11 Initial version
/// @date 2017-07-26 Rewrite with move semantics
/// @author Andreas Salzburger
/// @author Moritz Kiehn <msmk@cern.ch>

#ifndef ACTFW_WHITEBOARD_H
#define ACTFW_WHITEBOARD_H

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

/// A container to store arbitrary objects with ownership transfer.
///
/// This is an append-only container that takes ownership of the objects
/// added to it. Once an object has been added, it can only be read but not
/// be modified. Trying to replace an existing object is considered an error.
/// Its lifetime is bound to the liftime of the white board.
class WhiteBoard
{
public:
  WhiteBoard(std::unique_ptr<const Acts::Logger> logger
             = Acts::getDefaultLogger("WhiteBoard", Acts::Logging::INFO));

  // A WhiteBoard holds unique elements and can not be copied
  WhiteBoard(const WhiteBoard& other) = delete;
  WhiteBoard&
  operator=(const WhiteBoard&)
      = delete;

  /// Store an object on the white board and transfer ownership.
  ///
  /// @param name Identifier to store it under
  /// @param object Movable reference to the transferable object
  /// @returns ProcessCode::SUCCESS if the object was stored successfully
  template <typename T>
  ProcessCode
  add(const std::string& name, T&& object);

  /// Get access to a stored object.
  ///
  /// @param[in] name Identifier for the object
  /// @param[out] object A pointer to the object or nullptr on error
  /// @returns ProcessCode::SUCCESS if the object was found
  template <typename T>
  ProcessCode
  get(const std::string& name, const T*& object) const;

private:
  // type-erased value holder for move-constructible types
  struct IHolder
  {
    virtual ~IHolder() = default;
    virtual const std::type_info&
    type() const = 0;
  };
  template <typename T,
            typename
            = std::enable_if_t<std::is_nothrow_move_constructible<T>::value>>
  struct HolderT : public IHolder
  {
    T value;

    HolderT(T&& v) : value(std::move(v)) {}
    const std::type_info&
    type() const
    {
      return typeid(T);
    }
  };

  std::unique_ptr<const Acts::Logger> m_logger;
  std::map<std::string, std::unique_ptr<IHolder>> m_store;

  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
};

}  // namespace FW

inline FW::WhiteBoard::WhiteBoard(std::unique_ptr<const Acts::Logger> logger)
  : m_logger(std::move(logger))
{
}

template <typename T>
inline FW::ProcessCode
FW::WhiteBoard::add(const std::string& name, T&& object)
{
  if (0 < m_store.count(name)) {
    ACTS_FATAL("Object '" << name << "' already exists");
    return ProcessCode::ABORT;
  }
  m_store.emplace(name, std::make_unique<HolderT<T>>(std::forward<T>(object)));
  ACTS_VERBOSE("Added object '" << name << "'");
  return ProcessCode::SUCCESS;
}

template <typename T>
inline FW::ProcessCode
FW::WhiteBoard::get(const std::string& name, const T*& object) const
{
  auto it = m_store.find(name);
  if (it == m_store.end()) {
    object = nullptr;
    ACTS_FATAL("Object '" << name << "' does not exists");
    return ProcessCode::ABORT;
  }
  const IHolder* holder = it->second.get();
  if (typeid(T) != holder->type()) {
    object = nullptr;
    ACTS_FATAL("Type missmatch for object '" << name << "'");
    return ProcessCode::ABORT;
  }
  object = &(reinterpret_cast<const HolderT<T>*>(holder)->value);
  ACTS_VERBOSE("Retrieved object '" << name << "'");
  return ProcessCode::SUCCESS;
}

#endif  // ACTFW_WHITEBOARD_H
