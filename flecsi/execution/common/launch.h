/*
    @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
   /@@/////  /@@          @@////@@ @@////// /@@
   /@@       /@@  @@@@@  @@    // /@@       /@@
   /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
   /@@////   /@@/@@@@@@@/@@       ////////@@/@@
   /@@       /@@/@@//// //@@    @@       /@@/@@
   /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
   //       ///  //////   //////  ////////  //

   Copyright (c) 2016, Los Alamos National Security, LLC
   All rights reserved.
                                                                              */
#pragma once

/*! @file */

#if !defined(__FLECSI_PRIVATE__)
#error Do not include this file directly
#else
#include <flecsi/utils/debruijn.h>
#endif

#include <bitset>

namespace flecsi {

/*!
  Bitmasks for launch types.

  \note This enumeration is not scoped so that users can do things
        like:
        \code
        launch_t l(single | leaf);
        \endcode
 */

enum launch_mask_t : size_t {
  single = 1 << 0,
  index = 1 << 1,
  leaf = 1 << 2,
  inner = 1 << 3,
  idempotent = 1 << 4
}; // enum launch_mask_t

namespace execution {

// This will be used by the task_hash_t type to create hash keys for
// task registration. If you add more launch flags below, you will need
// to increase the launch_bits accordingly, i.e., launch_bits must
// be greater than or equal to the number of bits in the bitset for
// launch_t below.

constexpr size_t launch_bits = 5;

/*!
  Use a std::bitset to store launch information.

  @note This will most likely use 4 bytes of data for efficiency.
 */

using launch_t = std::bitset<launch_bits>;

/*!
  Enumeration of various task launch types. Not all of these may be
  supported by all runtimes. Unsupported launch information will be
  ignored.
 */

enum class launch_type_t : size_t {
  single,
  index,
  leaf,
  inner,
  idempotent
}; // enum launch_type_t

/*!
  Convert a processor mask to a processor type.
 */

inline launch_type_t
mask_to_type(launch_mask_t m) {
  return static_cast<launch_type_t>(flecsi::utils::debruijn32_t::index(m));
} // mask_to_type

/*!
  Macro to create repetitive interfaces.
 */

#define test_boolean_interface(name)                                           \
  inline bool launch_##name(const launch_t & l) {                              \
    return l.test(static_cast<size_t>(launch_type_t::name));                   \
  }

// clang-format off
test_boolean_interface(single)
test_boolean_interface(index)
test_boolean_interface(leaf)
test_boolean_interface(inner)
test_boolean_interface(idempotent)
// clang-format on

#undef test_boolean_interface

} // namespace execution
} // namespace flecsi
