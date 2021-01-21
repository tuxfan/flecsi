/*
   Copyright (c) 2016, Triad National Security, LLC
   All rights reserved.
                                                                              */
#pragma once

#include <flecsi/execution.hh>

#include "../data.hh"

namespace standalone {
namespace task {

using namespace flecsi;

void
advance_mpi(single<mpi_data *>::accessor<rw> data) {

  std::stringstream ss;
  ss << "mpi data: ";
  for(std::size_t i{0}; i < 10; ++i) {
    data->values[i] += 1;
    ss << data->values[i] << " ";
  } // for
  flog(info) << ss.str() << std::endl;
} // advance_mpi

} // namespace task
} // namespace standalone
