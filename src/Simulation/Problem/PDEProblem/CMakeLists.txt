add_subdirectory("DiscretizeInTime")
add_subdirectory("PDEAlgorithms")

SET(HEADERS_PDEProblem
  Simulation/Problem/PDEProblem/PDEProblem.hpp
  ${HEADERS_DiscretizeInTime}
  ${HEADERS_PDEAlgorithms}
CACHE INTERNAL "")

SET(SOURCES_PDEProblem
  Simulation/Problem/PDEProblem/PDEProblem.cpp
  ${SOURCES_DiscretizeInTime}
  ${SOURCES_PDEAlgorithms}
CACHE INTERNAL "")
