/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "specialization/control.hh"
#include "state.hh"
#include "tasks/norm.hh"

#include <flecsi/execution.hh>

#include <cmath>

namespace poisson {
namespace action {
using namespace flecsi;

int
analysis() {
  double sum = reduce<task::diff, exec::fold::sum>(m, ud(m), sd(m)).get();
  sum = execute<task::scale>(m, sum).get();
  const double l2 = sqrt(sum);
  flog(info) << "l2 error: " << l2 << std::endl;
  return 0;
} // analysis

control::action<analysis, cp::analysis> analysis_action;

} // namespace action
} // namespace poisson
