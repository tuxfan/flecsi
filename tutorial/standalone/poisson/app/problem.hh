/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "initialize.hh"
#include "options.hh"
#include "specialization/control.hh"
#include "state.hh"
#include "tasks/init.hh"

namespace poisson {
namespace action {
using namespace flecsi;

int
problem() {
  execute<task::eggcarton>(m, ud(m), fd(m), sd(m));
  return 0;
} // init_mesh

control::action<problem, cp::initialize> problem_action;
const auto problem_dep = problem_action.add(init_mesh_action);

} // namespace action
} // namespace poisson
