add_subdirectory("STKMesh")
add_subdirectory("Entity")
add_subdirectory("MeshQuery")
add_subdirectory("MeshPart")
add_subdirectory("MeshUtility")

SET(HEADERS_Mesh
  Mesh/Mesh.hpp
  ${HEADERS_STKMesh}
  ${HEADERS_Entity}
  ${HEADERS_MeshQuery}
  ${HEADERS_MeshPart}
  ${HEADERS_MeshUtility}
CACHE INTERNAL "")

SET(SOURCES_Mesh
  Mesh/Mesh.cpp
  ${SOURCES_STKMesh}
  ${SOURCES_Entity}
  ${SOURCES_MeshQuery}
  ${SOURCES_MeshPart}
  ${SOURCES_MeshUtility}
CACHE INTERNAL "")
