#include "Parser/InputParser/Datafile.hpp"
#include "Mesh/MeshPart/BoundaryMeshPart.hpp"
#include "Mesh/Entity/Entity.hpp"

namespace NSMesh {

BoundaryMeshPart::BoundaryMeshPart() {
  _entityInfoPerTopology.resize( 4 );
}


void
BoundaryMeshPart::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  _name = p_sectionName;
  
  // Read the id
  _id = p_datafile.Get<int>( p_sectionPath + "/" + p_sectionName + "/ID" );
}


void
BoundaryMeshPart::Initialize( int p_ndim ) {
  if( p_ndim == 2 ) {
    _primaryTopology = NSMesh::Topology::tp_EDGE;
  }
  else if( p_ndim == 3 ) {
    _primaryTopology = NSMesh::Topology::tp_FACE;
  }
  else {
    NSExceptions::InvalidArgument( true, "The boundary mesh part is not supported for the problem dimension used" );
  }
}
  

const int&
BoundaryMeshPart::GetId() const {
  return _id;
}


const std::string&
BoundaryMeshPart::GetName() const {
  return _name;
}


const NSMesh::Topology&
BoundaryMeshPart::GetPrimaryTopology() const {
  return _primaryTopology;
}


void
BoundaryMeshPart::SetElementTopologyInfo( int p_numVerticesPerElement, const NSMesh::ElementTopology& p_elementTopology ) {
  _numberOfVerticesPerElement = p_numVerticesPerElement;
  _elementTopology = p_elementTopology;
}


void
BoundaryMeshPart::AllocateMemoryPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, size_t p_size ) {
  // Get the entity container map for the specified topology
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  
  // Check that no memory has been allocated before
  NSExceptions::LogicError( entity_container_map.find( p_interiorPartName ) != entity_container_map.end(), "Memory has been allocated for the boundary mesh part" );
  
  // Allocate memory
  entity_container_map.insert( std::make_pair( p_interiorPartName, std::make_unique<NSMesh::NSEntity::EntityContainer<NSMesh::NSEntity::Entity*> >( p_topology ) ) );
  entity_container_map.find( p_interiorPartName )->second->AllocateSizeForEntity( p_size );
}


void
BoundaryMeshPart::SetEntityStatsPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName,
					     const NSMesh::Ownership& p_ownership, size_t p_data ) {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  entity_container_map.find( p_interiorPartName )->second->SetNumberOfEntity( p_ownership, p_data );
}


void
BoundaryMeshPart::SetGlobalNumEntityPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, size_t p_data ) {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  entity_container_map.find( p_interiorPartName )->second->SetGlobalNumberOfEntity( p_data );
}


void
BoundaryMeshPart::SetEntityOwnershipIndexesPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName,
							const NSMesh::Ownership& p_ownership, size_t p_size, size_t p_startIndex ) {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  entity_container_map.find( p_interiorPartName )->second->SetEntityOwnershipIndexes( p_ownership, p_size, p_startIndex );
}


void
BoundaryMeshPart::AddEntityAtIndex( NSMesh::NSEntity::Entity* p_entity, const std::string& p_interiorPartName,
				    const NSMesh::Topology& p_topology, unsigned p_gid, size_t p_index ) {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  entity_container_map.find( p_interiorPartName )->second->AddEntityAtIndex( p_entity, p_gid, p_index );
}


void
BoundaryMeshPart::RebuildGIDToLIDMapPerTopology( const NSMesh::Topology& p_topology ) {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  for( auto& entity_container_per_part : entity_container_map ) {
    entity_container_per_part.second->RebuildGIDToLIDMap();
  }
}


const std::vector<std::pair<long unsigned, long unsigned> >&
BoundaryMeshPart::GetEntityIndexVec( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, const NSMesh::Ownership& p_ownership ) const {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  return entity_container_map.find( p_interiorPartName )->second->GetEntityIndexVec( p_ownership );
}


const std::vector<NSMesh::NSEntity::Entity*>&
BoundaryMeshPart::GetEntityVec( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName ) const {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  return entity_container_map.find( p_interiorPartName )->second->GetEntityVec();
}


unsigned
BoundaryMeshPart::GetNumberOfEntity( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, const std::vector<NSMesh::Ownership>& p_ownership ) const {
  auto& entity_container_map = _entityInfoPerTopology[ p_topology ];
  return entity_container_map.find( p_interiorPartName )->second->GetNumberOfEntity( p_ownership );
}


NSMesh::ElementTopology
BoundaryMeshPart::GetMeshElementTopology() const {
  return _elementTopology;
}


int
BoundaryMeshPart::GetNumberOfVerticesPerElement() const {
  return _numberOfVerticesPerElement;
}


void
BoundaryMeshPart::SetAttachedInteriorMeshPartsName( const std::vector<std::string>& p_attachedInteriorMeshPartsName ) {
  _attachedInteriorMeshPartNames = p_attachedInteriorMeshPartsName;
}
  

const std::vector<std::string>&
BoundaryMeshPart::GetAttachedInteriorMeshPartsNames() const {
  return _attachedInteriorMeshPartNames;
}
  
} // namespace mesh
