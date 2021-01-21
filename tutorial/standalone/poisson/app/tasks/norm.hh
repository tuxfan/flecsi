/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "../../specialization/mesh.hh"

namespace poisson {
namespace task {
using namespace flecsi;

double
diff(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba) {
  auto a = m.mdspan<mesh::vertices>(aa);
  auto b = m.mdspan<mesh::vertices>(ba);

  double sum{0};
  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      sum += pow(a[j][i] - b[j][i], 2);
    } // for
  } // for

  return sum;
} // diff

double
scale(mesh::accessor<ro> m, double sum) {
  return pow(m.delta(), 2)*sum;
} // scale

void
discrete_operator(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<rw, ro> Aua) {
  auto u = m.mdspan<mesh::vertices>(ua);
  auto Au = m.mdspan<mesh::vertices>(Aua);

  const double w = 1.0 / pow(m.delta(), 2);

  // clang-format off
  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      Au[j][i] = w * (4.0 * u[j][i] -
        u[j][i + 1] - u[j][i - 1] - u[j + 1][i] - u[j - 1][i]);
    } // for
  } // for
  // clang-format on
} // residual

} // namespace task
} // namespace poisson
