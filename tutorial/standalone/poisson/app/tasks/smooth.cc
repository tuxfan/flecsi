/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/smooth.hh"

using namespace flecsi;

void
poisson::task::smooth(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa) {
  auto u = m.mdspan<mesh::vertices>(ua);
  auto f = m.mdspan<mesh::vertices>(fa);
  const auto dsqr = pow(m.delta(), 2);

  // clang-format off
  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      u[j][i] = 0.25 * (dsqr * f[j][i] +
        u[j][i + 1] + u[j][i - 1] + u[j + 1][i] + u[j - 1][i]);
    } // for
  } // for
  // clang format on
} // smooth

void poisson::task::jacobi(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua_new,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<ro, ro> fa,
  double omega) {
  auto u_new = m.mdspan<mesh::vertices>(ua_new);
  auto u = m.mdspan<mesh::vertices>(ua);
  auto f = m.mdspan<mesh::vertices>(fa);
  const auto dsqr = pow(m.delta(), 2);

  // clang-format off
  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      u_new[j][i] = (1.0 - omega) * u[j][i] + omega*(0.25 * (dsqr * f[j][i] +
        u[j][i + 1] + u[j][i - 1] + u[j + 1][i] + u[j - 1][i]));
    } // for
  } // for
  // clang format on
} // jacobi
