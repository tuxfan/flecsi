#------------------------------------------------------------------------------#
# Copyright (c) 2016, Triad National Security, LLC
# All rights reserved
#------------------------------------------------------------------------------#

include_directories(${PROJECT_SOURCE_DIR})

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#

add_executable(poisson
  poisson.cc
  analysis.hh
  initialize.hh
  options.hh
  problem.hh
  solve.hh
  state.hh
  tasks/init.hh
  tasks/norm.hh
  tasks/smooth.hh
  ${PROJECT_SOURCE_DIR}/specialization/control.hh
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#

target_link_libraries(poisson
  ${FleCSI_LIBRARIES}
  ${FleCSI_LIBRARY_DEPENDENCIES}
)
