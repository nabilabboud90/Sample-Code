add_subdirectory("Problem")

SET(HEADERS_Simulation
  Simulation/Simulation.hpp
  ${HEADERS_Problem}
CACHE INTERNAL "")

SET(SOURCES_Simulation
  Simulation/Simulation.cpp
  ${SOURCES_Problem}
CACHE INTERNAL "")
