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

#include "flecsi/data/accessor.hh"
#include "flecsi/data/copy_plan.hh"
#include "flecsi/data/layout.hh"
#include "flecsi/data/privilege.hh"
#include "flecsi/data/topology.hh"
#include "flecsi/flog.hh"
#include "flecsi/topo/core.hh"
#include "flecsi/topo/index.hh"
#include "flecsi/topo/narray/types.hh"
#include "flecsi/topo/utility_types.hh"
#include "flecsi/util/array_ref.hh"

#include <utility>

namespace flecsi {
namespace topo {

/*----------------------------------------------------------------------------*
  Narray Topology.
 *----------------------------------------------------------------------------*/

template<typename Policy>
struct narray : narray_base, with_ragged<Policy>, with_meta<Policy> {
  using index_space = typename Policy::index_space;
  using index_spaces = typename Policy::index_spaces;
  using axis = typename Policy::axis;
  using axes = typename Policy::axes;
  using coord = narray_impl::coord;
  using hypercube = narray_impl::hypercube;
  using coloring_definition = narray_impl::coloring_definition;
  using id = util::id;

  static constexpr std::size_t dimension = Policy::dimension;

  template<std::size_t>
  struct access;

  narray(coloring const & c)
    : with_ragged<Policy>(c.colors), with_meta<Policy>(c.colors),
      part_(make_partitions(c,
        index_spaces(),
        std::make_index_sequence<index_spaces::size>())),
      plan_(make_plans(c,
        index_spaces(),
        std::make_index_sequence<index_spaces::size>())),
      meta_(c.colors) {
    init_ragged(index_spaces());
    init_meta(c);
  }

  struct meta_data {
    using scoord = std::array<std::size_t, dimension>;
    using shypercube = std::array<scoord, 2>;
    std::uint32_t faces;

    std::array<std::array<std::size_t, dimension>, index_spaces::size> global;
    std::array<std::array<std::size_t, dimension>, index_spaces::size> offset;
    std::array<std::array<std::size_t, dimension>, index_spaces::size> extents;
    std::array<shypercube, index_spaces::size> logical;
    std::array<shypercube, index_spaces::size> extended;

    static_assert(std::is_trivial<typename Policy::meta_data>::value,
      "illegal meta data type: must be a trivial type");
    typename Policy::meta_data meta;
  };

  static inline const typename field<meta_data,
    data::single>::template definition<meta<narray<Policy>>>
    meta_field;

  util::key_array<repartitioned, index_spaces> part_;
  util::key_array<data::copy_plan, index_spaces> plan_;
  typename meta<narray<Policy>>::core meta_;

  std::size_t colors() const {
    return part_.front().colors();
  }

  template<index_space S>
  data::region & get_region() {
    return part_.template get<S>();
  }

  template<index_space S>
  const data::partition & get_partition(field_id_t) const {
    return part_.template get<S>();
  }

  template<typename Type,
    data::layout Layout,
    typename Topo,
    typename Topo::index_space Space>
  void ghost_copy(data::field_reference<Type, Layout, Topo, Space> const & f) {
    plan_.template get<Space>().issue_copy(f.fid());
  }

private:
  template<auto... Value, std::size_t... Index>
  util::key_array<repartitioned, util::constants<Value...>> make_partitions(
    narray_base::coloring const & c,
    util::constants<Value...> /* index spaces to deduce pack */,
    std::index_sequence<Index...>) {
    flog_assert(c.idx_colorings.size() == sizeof...(Value),
      c.idx_colorings.size()
        << " sizes for " << sizeof...(Value) << " index spaces");
    return {{make_repartitioned<Policy, Value>(
      c.colors, make_partial<idx_size>(c.idx_colorings[Index]))...}};
  }

  template<index_space S>
  data::copy_plan make_plan(index_coloring const & ic) {
    std::vector<std::size_t> num_intervals;

    execute<idx_itvls, mpi>(ic, num_intervals);

    // clang-format off
    auto dest_task = [&ic](auto f) {
      execute<set_dests, mpi>(f, ic.intervals);
    };

    auto ptrs_task = [&ic](auto f) {
      execute<set_ptrs<Policy::template privilege_count<S>>, mpi>(f, ic.points);
    };

    return {*this, num_intervals, dest_task, ptrs_task, util::constant<S>()};
    // clang-format on
  }

  template<auto... Value, std::size_t... Index>
  util::key_array<data::copy_plan, util::constants<Value...>> make_plans(
    narray_base::coloring const & c,
    util::constants<Value...> /* index spaces to deduce pack */,
    std::index_sequence<Index...>) {
    flog_assert(c.idx_colorings.size() == sizeof...(Value),
      c.idx_colorings.size()
        << " sizes for " << sizeof...(Value) << " index spaces");
    return {{make_plan<Value>(c.idx_colorings[Index])...}};
  }

  static void set_meta(
    typename field<meta_data, data::single>::template accessor<wo> m,
    narray_base::coloring const & c) {
    meta_data & md = m;

    for(std::size_t i{0}; i < index_spaces::size; ++i) {
      auto const & cg = c.idx_colorings[i].global;
      auto const & co = c.idx_colorings[i].offset;
      auto const & ce = c.idx_colorings[i].extents;

      flog_assert(cg.size() == dimension,
        "invalid #axes(" << cg.size() << ") must be: " << dimension);
      flog_assert(co.size() == dimension,
        "invalid #axes(" << co.size() << ") must be: " << dimension);
      flog_assert(ce.size() == dimension,
        "invalid #axes(" << ce.size() << ") must be: " << dimension);

      md.faces = c.idx_colorings[i].faces;

      std::copy_n(cg.begin(), dimension, md.global[i].begin());
      std::copy_n(co.begin(), dimension, md.offset[i].begin());
      std::copy_n(ce.begin(), dimension, md.extents[i].begin());

      auto const & cl = c.idx_colorings[i].logical;
      std::copy_n(cl[0].begin(), dimension, md.logical[i][0].begin());
      std::copy_n(cl[1].begin(), dimension, md.logical[i][1].begin());

      auto const & cex = c.idx_colorings[i].extended;
      std::copy_n(cex[0].begin(), dimension, md.extended[i][0].begin());
      std::copy_n(cex[1].begin(), dimension, md.extended[i][1].begin());
    } // for
  } // set_meta

  void init_meta(narray_base::coloring const & c) {
    execute<set_meta, mpi>(meta_field(this->meta_), c);
  }

  template<index_space... SS>
  void init_ragged(util::constants<SS...>) {
    (this->template extend_offsets<SS>(), ...);
  }
}; // struct narray

/*----------------------------------------------------------------------------*
  Narray Access.
 *----------------------------------------------------------------------------*/

template<typename Policy>
template<std::size_t Privileges>
struct narray<Policy>::access {
  template<const auto & F>
  using accessor = data::accessor_member<F, Privileges>;
  util::key_array<resize::accessor<ro>, index_spaces> size_;

  accessor<meta_field> meta_;

  access() {}

  enum class range : std::size_t {
    logical,
    extended,
    all,
    boundary_low,
    boundary_high,
    ghost_low,
    ghost_high,
    global
  };

  using hypercubes = index::has<range::logical,
    range::extended,
    range::all,
    range::boundary_low,
    range::boundary_high,
    range::ghost_low,
    range::ghost_high,
    range::global>;

  template<axis A>
  bool is_low() {
    return (meta_.get().faces >> A * 2) & narray_impl::low;
  }

  template<axis A>
  bool is_high() {
    return (meta_.get().faces >> A * 2) & narray_impl::high;
  }

  template<axis A>
  bool is_interior() {
    return !is_low<A>() && !is_high<A>();
  }

  // This should be protected
  auto & meta() {
    return meta_.get().meta;
  }

  template<index_space S, axis A, range SE>
  std::size_t size() {
    auto const & md = meta_.get();
    static_assert(
      std::size_t(SE) < hypercubes::size, "invalid size identifier");
    if constexpr(SE == range::logical) {
      return md.logical[S][1][A];
    }
    else if constexpr(SE == range::extended) {
      return md.extended[S][1][A];
    }
    else if constexpr(SE == range::all) {
      return md.extents[S][A];
    }
    else if constexpr(SE == range::boundary_low) {
      return md.logical[S][0][A] - md.extended[S][0][A];
    }
    else if constexpr(SE == range::boundary_high) {
      return md.extended[S][1][A] - md.logical[S][1][A];
    }
    else if constexpr(SE == range::ghost_low) {
      return md.logical[S][0][A];
    }
    else if constexpr(SE == range::ghost_high) {
      return md.extents[S][A] - md.logical[S][1][A];
    }
    else if constexpr(SE == range::global) {
      return md.global[S][A];
    }
  }

  template<index_space S, axis A, range SE>
  auto extents() {
    auto const & md = meta_.get();
    static_assert(
      std::size_t(SE) < hypercubes::size, "invalid extents identifier");
    if constexpr(SE == range::logical) {
      return make_ids<S>(
        util::iota_view<util::id>(md.logical[S][0][A], md.logical[S][1][A]));
    }
    else if constexpr(SE == range::extended) {
      return make_ids<S>(
        util::iota_view<util::id>(md.extended[S][0][A], md.extended[S][1][A]));
    }
    else if constexpr(SE == range::all) {
      return make_ids<S>(util::iota_view<util::id>(0, md.extents[S][A]));
    }
    else if constexpr(SE == range::boundary_low) {
      return make_ids<S>(util::iota_view<util::id>(
        md.extended[S][0][A], md.extended[S][0][A] + size<S, A, SE>()));
    }
    else if constexpr(SE == range::boundary_high) {
      return make_ids<S>(util::iota_view<util::id>(
        md.extended[S][1][A], md.extended[S][1][A] + size<S, A, SE>()));
    }
    else if constexpr(SE == range::ghost_low) {
      return make_ids<S>(util::iota_view<util::id>(
        md.logical[S][0][A], md.logical[S][0][A] + size<S, A, SE>()));
    }
    else if constexpr(SE == range::ghost_high) {
      return make_ids<S>(util::iota_view<util::id>(
        md.logical[S][1][A], md.logical[S][1][A] + size<S, A, SE>()));
    }
    else {
      flog_error("invalid range");
    }
  }

  template<index_space S, axis A, range SE>
  std::size_t offset() {
    auto const & md = meta_.get();
    static_assert(
      std::size_t(SE) < hypercubes::size, "invalid offset identifier");
    if constexpr(SE == range::logical) {
      return md.logical[S][0][A];
    }
    else if constexpr(SE == range::extended) {
      return md.extended[S][0][A];
    }
    else if constexpr(SE == range::all) {
      return md.extents[S][A];
    }
    else if constexpr(SE == range::boundary_low) {
      return md.extended[S][0][A];
    }
    else if constexpr(SE == range::boundary_high) {
      return md.logical[S][1][A];
    }
    else if constexpr(SE == range::ghost_low) {
      return 0;
    }
    else if constexpr(SE == range::ghost_high) {
      return md.logical[S][1][A];
    }
    else if constexpr(SE == range::global) {
      return md.offset[S][A];
    }
  }

  template<index_space S, typename T, std::size_t P>
  auto mdspan(data::accessor<data::dense, T, P> const & a) {
    auto const s = a.span();
    return util::mdspan<typename decltype(s)::element_type, dimension>(
      s.data(), meta_.get().extents[S]);
  }

  template<class F>
  void send(F && f) {
    std::size_t i{0};
    for(auto & a : size_) {
      f(a, [&i](typename Policy::slot & n) { return n->part_[i++].sizes(); });
    }
    meta_.topology_send(f, &narray::meta_);
  }
}; // struct narray<Policy>::access

/*----------------------------------------------------------------------------*
  Define Base.
 *----------------------------------------------------------------------------*/

template<>
struct detail::base<narray> {
  using type = narray_base;
}; // struct detail::base<narray>

} // namespace topo
} // namespace flecsi