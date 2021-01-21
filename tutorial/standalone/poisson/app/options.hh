/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "flecsi/execution.hh"

flecsi::program_option<std::size_t>
  x_extents("x-extents", "The x extents of the mesh.", 1);
flecsi::program_option<std::size_t>
  y_extents("y-extents", "The y extents of the mesh.", 1);
