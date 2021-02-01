/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#include "analysis.hh"
#include "state.hh"
#include "tasks/norm.hh"

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

#include <cmath>

using namespace flecsi;

int
poisson::action::analysis() {
  double sum = reduce<task::diff, exec::fold::sum>(m, ud(m), sd(m)).get();
  sum = execute<task::scale>(m, sum).get();
  const double l2 = sqrt(sum);
  flog(info) << "l2 error: " << l2 << std::endl;
  return 0;
} // analysis
