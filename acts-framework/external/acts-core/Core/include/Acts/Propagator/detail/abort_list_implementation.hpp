// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include "Acts/Propagator/detail/condition_uses_result_type.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {
    /// This is the caller used if the condition uses the result
    /// and that result type exists
    template <bool has_result = true>
    struct condition_caller
    {
      template <typename condition,
                typename result_t,
                typename propagator_state_t,
                typename stepper_t>
      static bool
      check(const condition&    c,
            const result_t&     r,
            propagator_state_t& state,
            const stepper_t&    stepper)
      {
        using action_type = action_type_t<condition>;
        using result_type = result_type_t<action_type>;

        return c(r.template get<result_type>(), state, stepper);
      }
    };

    /// This is the caller used if the condition only uses the cache
    /// it has no access to the result
    template <>
    struct condition_caller<false>
    {
      template <typename condition,
                typename result_t,
                typename propagator_state_t,
                typename stepper_t>
      static bool
      check(const condition& c,
            const result_t& /*result*/,
            propagator_state_t& state,
            const stepper_t&    stepper)
      {
        return c(state, stepper);
      }
    };
  }  // end of anonymous namespace

  template <typename... conditions>
  struct abort_list_impl;

  /// This is the check call on the a list of conditions
  /// it calls the aparant condition and broadcasts
  /// the call to the remaining ones
  template <typename first, typename... others>
  struct abort_list_impl<first, others...>
  {
    template <typename T,
              typename result_t,
              typename propagator_state_t,
              typename stepper_t>
    static bool
    check(const T&            conditions_tuple,
          const result_t&     result,
          propagator_state_t& state,
          const stepper_t&    stepper)
    {

      // get the right helper for calling the abort condition
      constexpr bool has_result = condition_uses_result_type<first>::value;
      using caller_type         = condition_caller<has_result>;

      // get the cache abort condition
      const auto& this_condition = std::get<first>(conditions_tuple);

      // - check abort conditions recursively
      // - make use of short-circuit evaluation
      // -> skip remaining conditions if this abort condition evaluates to true
      bool abort = caller_type::check(this_condition, result, state, stepper)
          || abort_list_impl<others...>::check(
                       conditions_tuple, result, state, stepper);

      return abort;
    }
  };

  /// This is the check call on the a last of all conditions
  template <typename last>
  struct abort_list_impl<last>
  {
    template <typename T,
              typename result_t,
              typename propagator_state_t,
              typename stepper_t>
    static bool
    check(const T&            conditions_tuple,
          const result_t&     result,
          propagator_state_t& state,
          const stepper_t&    stepper)
    {
      // get the right helper for calling the abort condition
      constexpr bool has_result     = condition_uses_result_type<last>::value;
      const auto&    this_condition = std::get<last>(conditions_tuple);

      return condition_caller<has_result>::check(
          this_condition, result, state, stepper);
    }
  };

  /// This is the empty call list - never abort
  template <>
  struct abort_list_impl<>
  {
    template <typename T,
              typename result_t,
              typename propagator_state_t,
              typename stepper_t>
    static bool
    check(const T& /*unused*/,
          const result_t& /*result*/,
          propagator_state_t& /*state*/,
          const stepper_t& /*unused*/)
    {
      return false;
    }
  };

}  // namespace detail

}  // namespace Acts
