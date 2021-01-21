/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "state.hh"
#include "tasks/io.hh"

#include <flecsi/flog.hh>

namespace poisson {
namespace action {

int
finalize() {
  execute<task::io, mpi>(m, sd(m));
  return 0;
} // init_mesh

control::action<finalize, cp::finalize> finalize_action;

} // namespace action
} // namespace poisson
