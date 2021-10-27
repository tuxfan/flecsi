/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/smooth.hh"

using namespace flecsi;

void
poisson::task::smooth(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa,
  bool red) {
  auto u = m.mdspan<mesh::vertices>(ua);
  auto f = m.mdspan<mesh::vertices>(fa);
  const auto dxdy = m.dxdy();
  const auto dx_over_dy = m.xdelta() / m.ydelta();
  const auto dy_over_dx = m.ydelta() / m.xdelta();
  const auto factor = 1.0 / (2 * (dx_over_dy + dy_over_dx));

  // clang-format off
  for(auto j : m.vertices<mesh::y_axis>()) {
    forall(i, red ? m.red(j) : m.black(j), "smooth") {
      u[j][i] = factor *
                (dxdy * f[j][i] +
                 dy_over_dx * (u[j][i + 1] + u[j][i - 1]) +
                 dx_over_dy * (u[j + 1][i] + u[j - 1][i]));
    }; // for
  } // for
  // clang format on
} // smooth
