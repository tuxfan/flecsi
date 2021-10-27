/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "tasks/init.hh"

#include <cmath>

using namespace flecsi;

constexpr double PI = 3.14159;
constexpr double K = 12.0;
constexpr double L = 2.0;

void
poisson::task::eggcarton(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  field<double>::accessor<wo, na> fa,
  field<double>::accessor<wo, na> sa,
  field<double>::accessor<wo, na> Aua) {
  auto u = m.mdspan<mesh::vertices>(ua);
  auto f = m.mdspan<mesh::vertices>(fa);
  auto s = m.mdspan<mesh::vertices>(sa);
  auto Au = m.mdspan<mesh::vertices>(Aua);
  const double sq_klpi = pow(PI, 2) * (pow(K, 2) + pow(L, 2));

  flog(info) << "dxdy: " << m.dxdy() << std::endl;

  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      f[j][i] = sq_klpi * sin(K * PI * x) * sin(L * PI * y);
      const double solution = sin(K * PI * x) * sin(L * PI * y);
      s[j][i] = solution;
      Au[j][i] = 0.0;
      u[j][i] = 0.0;
    } // for
  } // for
} // eggcarton

void
poisson::task::constant(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double value) {
  auto f = m.mdspan<mesh::vertices>(fa);
  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      f[j][i] = value;
    } // for
  } // for
}

void
poisson::task::redblack(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa) {
  auto f = m.mdspan<mesh::vertices>(fa);
  for(auto j : m.vertices<mesh::y_axis, mesh::interior>()) {
    forall(i, m.red(j), "red") {
      f[j][i] = m.global_id<mesh::x_axis>(i);
    };
    forall(i, m.black(j), "black") {
      f[j][i] = -int(m.global_id<mesh::x_axis>(i));
    };
  } // for
}
