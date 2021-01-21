/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "options.hh"
#include "specialization/mesh.hh"

namespace poisson {

const flecsi::field<double>::definition<mesh, mesh::vertices> ud;
const flecsi::field<double>::definition<mesh, mesh::vertices> fd;
const flecsi::field<double>::definition<mesh, mesh::vertices> sd;
const flecsi::field<double>::definition<mesh, mesh::vertices> Aud;

mesh::slot m;
mesh::cslot coloring;

} // namespace poisson
