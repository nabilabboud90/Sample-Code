add_subdirectory("PDEProblem")

SET(HEADERS_Problem
  Simulation/Problem/Problem.hpp
  Simulation/Problem/Algorithm.hpp
  ${HEADERS_PDEProblem}
CACHE INTERNAL "")

SET(SOURCES_Problem
  Simulation/Problem/Problem.cpp
  Simulation/Problem/Algorithm.cpp
  ${SOURCES_PDEProblem}
CACHE INTERNAL "")
