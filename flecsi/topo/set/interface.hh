/*
    @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
   /@@/////  /@@          @@////@@ @@////// /@@
   /@@       /@@  @@@@@  @@    // /@@       /@@
   /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
   /@@////   /@@/@@@@@@@/@@       ////////@@/@@
   /@@       /@@/@@//// //@@    @@       /@@/@@
   /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
   //       ///  //////   //////  ////////  //

   Copyright (c) 2016, Triad National Security, LLC
   All rights reserved.
                                                                              */
#pragma once

/*! @file */

#if !defined(__FLECSI_PRIVATE__)
#error Do not include this file directly!
#endif
#include "flecsi/data/reference.hh"
#include "flecsi/topo/core.hh" // base
#include "flecsi/topo/set/types.hh"

namespace flecsi {
namespace topo {

template<typename Policy>
struct set : set_base, with_meta<Policy> {

  template<std::size_t>
  struct access;

  set(coloring const & c)
    : with_meta<Policy>(c.colors) {
  }

}; // struct set

template<typename Policy>
template<std::size_t Privileges>
struct set<Policy>::access {
}; // struct set<Policy>::access

template<>
struct detail::base<set> {
  using type = set_base;
};

} // namespace topo
} // namespace flecsi
