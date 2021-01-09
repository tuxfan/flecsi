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

#include "flecsi/data/topology.hh"
#include "flecsi/execution.hh"
#include "flecsi/topo/index.hh"
#include "flecsi/util/color_map.hh"
#include "flecsi/util/mpi.hh"
#include "flecsi/util/serialize.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <vector>

namespace flecsi {
namespace topo {
namespace narray_impl {

enum masks : uint32_t { interior = 0b00, low = 0b01, high = 0b10 };

enum axes : std::size_t { x_axis, y_axis, z_axis };

using coord = std::vector<std::size_t>;
// using rect = std::vector<coord>; // can't be std::array<coord, 2>
using rect = std::array<coord, 2>;
using interval = std::pair<std::size_t, std::size_t>;

struct index_coloring {

  /*
    The local extents of this color.
   */

  coord extents;

  /*
    The global coordinate offset of the local block.
    Local to global id translation can be computed with this.
   */

  coord global;

  /*
    Owned always span a single subregion.
   */

  rect owned;

  /*
    Exclusive always span a single rect.
   */

  rect exclusive;

  /*
    Boolean indicating whether or not a copy plan should be created.
   */

  bool create_plan = true;

  /*
    Offsets on the remote color.
   */

  std::map<std::size_t, /* over colors */
    std::vector<std::pair</* local ghost offset, remote shared offset */
      std::size_t,
      std::size_t>>>
    points;

  /*
    Local ghost intervals.
   */

  std::vector<std::pair<std::size_t, std::size_t>> intervals;
};

} // namespace narray_impl

/*----------------------------------------------------------------------------*
  Base.
 *----------------------------------------------------------------------------*/

struct narray_base {
  using index_coloring = narray_impl::index_coloring;

  struct coloring {
    std::size_t colors;
    std::vector<index_coloring> idx_colorings;
  }; // struct coloring

  static std::size_t idx_size(index_coloring const & ic, std::size_t) {
    std::size_t allocation{1};
    for(auto e : ic.extents) {
      allocation *= e;
    }
    return allocation;
  }

  static void idx_itvls(index_coloring const & ic,
    std::vector<std::size_t> & num_intervals) {
    num_intervals =
      util::mpi::all_gather([&ic](int, int) { return ic.intervals.size(); });
  }

  static void set_dests(field<data::intervals::Value>::accessor<wo> a,
    std::vector<std::pair<std::size_t, std::size_t>> const & intervals) {
    flog_assert(a.span().size() == intervals.size(), "interval size mismatch");
    std::size_t i{0};
    for(auto it : intervals) {
#if 0
      flog(warn) << "a[" << i << "] = intervals::make({" << it.first << ", "
                 << it.second << "}, " << process() << ")" << std::endl;
#endif
      a[i++] = data::intervals::make({it.first, it.second}, process());
    } // for
  }

  static void set_ptrs(field<data::points::Value>::accessor<wo> a,
    std::map<std::size_t,
      std::vector<std::pair<std::size_t, std::size_t>>> const & shared_ptrs) {
    for(auto const & si : shared_ptrs) {
      for(auto p : si.second) {
        // si.first: owner
        // p.first: local ghost offset
        // p.second: remote shared offset
        a[p.first] = data::points::make(si.first, p.second);
#if 0
        flog(warn) << "a[" << p.first << "] = points::make(" << si.first << ", "
                   << p.second << ")" << std::endl;
#endif
      } // for
    } // for
  }
}; // struct narray_base

} // namespace topo

/*----------------------------------------------------------------------------*
  Serialization Rules
 *----------------------------------------------------------------------------*/

template<>
struct util::serial<topo::narray_impl::index_coloring> {
  using type = topo::narray_impl::index_coloring;
  template<class P>
  static void put(P & p, const type & s) {
    serial_put(p,
      std::tie(s.extents,
        s.global,
        s.owned,
        s.exclusive,
        s.create_plan,
        s.points,
        s.intervals));
  }
  static type get(const std::byte *& p) {
    const serial_cast r{p};
    return type{r, r, r, r, r, r, r};
  }
};

} // namespace flecsi
