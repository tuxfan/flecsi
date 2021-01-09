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

#define __FLECSI_PRIVATE__
#include "flecsi/data.hh"
#include "flecsi/execution.hh"
#include "flecsi/flog.hh"
#include "flecsi/topo/narray/coloring_utils.hh"
#include "flecsi/topo/narray/interface.hh"
#include "flecsi/util/unit.hh"

using namespace flecsi;

struct mesh : topo::specialization<topo::narray, mesh> {

  enum index_space { entities };
  using index_spaces = has<entities>;
  enum rect { owned, exclusive, all };
  enum axis { x_axis, y_axis };
  using axes = has<x_axis, y_axis>;

  using coord = topo::narray_impl::coord;
  using coloring_definition = topo::narray_impl::coloring_definition;

  static constexpr std::size_t dimension = 2;

  template<auto>
  static constexpr std::size_t privilege_count = 2;

  /*--------------------------------------------------------------------------*
    Interface
   *--------------------------------------------------------------------------*/

  template<class B>
  struct interface : B {

    template<axis A, rect R = owned>
    std::size_t size() {
      switch(R) {
        case owned:
          return B::template size<index_space::entities, A, B::rect::owned>();
          break;
        case exclusive:
          return B::
            template size<index_space::entities, A, B::rect::exclusive>();
          break;
        case all:
          return B::template size<index_space::entities, A, B::rect::all>();
          break;
      }
    }

    template<axis A, rect R = owned>
    auto extents() {
      switch(R) {
        case owned:
          return B::
            template extents<index_space::entities, A, B::rect::owned>();
          break;
        case exclusive:
          return B::
            template extents<index_space::entities, A, B::rect::exclusive>();
          break;
        case all:
          return B::template extents<index_space::entities, A, B::rect::all>();
          break;
      }
    }
  };

  static coloring color(std::vector<coloring_definition> index_definitions) {
    auto [colors, index_colorings] =
      topo::narray_impl::color<2, 1>(index_definitions);

    flog_assert(colors == processes(),
      "current implementation is restricted to 1-to-1 mapping");

    coloring c;
    c.colors = colors;
    for(auto idx : index_colorings) {
      for(auto ic : idx) {
        c.idx_colorings.emplace_back(ic.second);
      }
    }
    return c;
  } // color
};

void
extents(mesh::accessor<ro> m, field<std::size_t>::accessor<wo, na> ca) {
  auto c = m.mdspan<mesh::entities>(ca);
  for(auto j : m.extents<mesh::y_axis>()) {
    for(auto i : m.extents<mesh::x_axis>()) {
      c[j][i] = color();
    } // for
  } // for
}

void
print(mesh::accessor<ro> m, field<std::size_t>::accessor<ro, ro> ca) {
  auto c = m.mdspan<mesh::entities>(ca);
  std::stringstream ss;
  for(int j{int(m.size<mesh::y_axis, mesh::all>() - 1)}; j >= 0; --j) {
    for(auto i : m.extents<mesh::x_axis, mesh::all>()) {
      ss << c[j][i] << " ";
    } // for
    ss << std::endl;
  } // for
  ss << std::endl;
  for(int j{int(m.size<mesh::y_axis, mesh::all>() - 1)}; j >= 0; --j) {
    for(auto i : m.extents<mesh::x_axis, mesh::all>()) {
      std::size_t id = i + j * int(m.size<mesh::x_axis, mesh::all>());
      ss << id << " ";
    } // for
    ss << std::endl;
  } // for
  flog(warn) << ss.str() << std::endl;
}

mesh::slot m;
mesh::cslot coloring;

const field<std::size_t>::definition<mesh, mesh::entities> cs;

int
narray_driver() {
  UNIT {
    // mesh::coord indices{8, 8};
    mesh::coord indices{25, 10};
    auto colors = topo::narray_impl::distribute(processes(), indices);
    flog(warn) << log::insert(colors, "colors") << std::endl;
    // mesh::coord colors{2, 2};
#if 1
    std::vector<mesh::coloring_definition> index_definitions = {
      {colors, indices}};
    coloring.allocate(index_definitions);
    m.allocate(coloring.get());
    execute<extents>(m, cs(m));
    execute<print>(m, cs(m));
#endif
  };
} // coloring_driver

flecsi::unit::driver<narray_driver> driver;
