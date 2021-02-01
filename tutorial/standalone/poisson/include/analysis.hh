/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "specialization/control.hh"

namespace poisson {
namespace action {

int analysis();
inline control::action<analysis, cp::analysis> analysis_action;

} // namespace action
} // namespace poisson
