#include "Parser/InputParser/Datafile.hpp"
#include "Mesh/MeshPart/InteriorMeshPart.hpp"
#include "Mesh/Entity/Entity.hpp"

namespace NSMesh {

InteriorMeshPart::InteriorMeshPart() {
  _entityInfoPerTopology.push_back( std::make_unique<NSMesh::NSEntity::EntityContainer<NSMesh::NSEntity::Entity*> >( NSMesh::Topology::tp_ELEMENT ) );
  _entityInfoPerTopology.push_back( std::make_unique<NSMesh::NSEntity::EntityContainer<NSMesh::NSEntity::Entity*> >( NSMesh::Topology::tp_FACE ) );
  _entityInfoPerTopology.push_back( std::make_unique<NSMesh::NSEntity::EntityContainer<NSMesh::NSEntity::Entity*> >( NSMesh::Topology::tp_EDGE ) );
  _entityInfoPerTopology.push_back( std::make_unique<NSMesh::NSEntity::EntityContainer<NSMesh::NSEntity::Entity*> >( NSMesh::Topology::tp_VERTEX ) );
}


void
InteriorMeshPart::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  _name = p_sectionName;
  
  // Read the id
  _id = p_datafile.Get<int>( p_sectionPath + "/" + p_sectionName + "/ID" );
}


const int&
InteriorMeshPart::GetId() const {
  return _id;
}


const std::string&
InteriorMeshPart::GetName() const {
  return _name;
}


NSMesh::Topology
InteriorMeshPart::GetPrimaryTopology() const {
  return NSMesh::Topology::tp_ELEMENT;
}


void
InteriorMeshPart::AllocateMemoryPerTopology( const NSMesh::Topology& p_topology, size_t p_size ) {
  _entityInfoPerTopology[ p_topology ]->AllocateSizeForEntity( p_size );
}


void
InteriorMeshPart::SetEntityStatsPerTopology( const NSMesh::Topology& p_topology, const NSMesh::Ownership& p_ownership, size_t p_data ) {
  _entityInfoPerTopology[ p_topology ]->SetNumberOfEntity( p_ownership, p_data );
}


void
InteriorMeshPart::SetGlobalNumEntityPerTopology( const NSMesh::Topology& p_topology, size_t p_data ) {
  _entityInfoPerTopology[ p_topology ]->SetGlobalNumberOfEntity( p_data );
}

  
void
InteriorMeshPart::SetEntityOwnershipIndexesPerTopology( const NSMesh::Topology& p_topology, const NSMesh::Ownership& p_ownership, size_t p_size, size_t p_startIndex ) {
  _entityInfoPerTopology[ p_topology ]->SetEntityOwnershipIndexes( p_ownership, p_size, p_startIndex );
}


void
InteriorMeshPart::AddEntityAtIndex( NSMesh::NSEntity::Entity* p_entity, const NSMesh::Topology& p_topology, unsigned p_gid, size_t p_index ) {
  _entityInfoPerTopology[ p_topology ]->AddEntityAtIndex( p_entity, p_gid, p_index );
}


void
InteriorMeshPart::RebuildGIDToLIDMapPerTopology( const NSMesh::Topology& p_topology ) {
  _entityInfoPerTopology[ p_topology ]->RebuildGIDToLIDMap();
}


const std::vector<std::pair<long unsigned, long unsigned> >&
InteriorMeshPart::GetEntityIndexVec( const NSMesh::Topology& p_topology, const NSMesh::Ownership& p_ownership ) const {
  return _entityInfoPerTopology[ p_topology ]->GetEntityIndexVec( p_ownership );
}


const std::vector<NSMesh::NSEntity::Entity*>&
InteriorMeshPart::GetEntityVec( const NSMesh::Topology& p_topology ) const {
  return _entityInfoPerTopology[ p_topology ]->GetEntityVec();
}


unsigned
InteriorMeshPart::GetNumberOfEntity( const NSMesh::Topology& p_topology, const std::vector<NSMesh::Ownership>& p_ownership ) const {
  return _entityInfoPerTopology[ p_topology ]->GetNumberOfEntity( p_ownership );
}


NSMesh::ElementTopology
InteriorMeshPart::GetMeshElementTopology() const {
  return _elementTopology;
}


int
InteriorMeshPart::GetNumberOfVerticesPerElement() const {
  return _numberOfVerticesPerElement;
}


void
InteriorMeshPart::SetElementTopologyInfo( int p_numVerticesPerElement, const NSMesh::ElementTopology& p_elementTopology ) {
  _numberOfVerticesPerElement = p_numVerticesPerElement;
  _elementTopology = p_elementTopology;
}


bool
InteriorMeshPart::IsPresent( unsigned p_gid, const NSMesh::Topology& p_topology ) const {
  return _entityInfoPerTopology[ p_topology ]->IsPresent( p_gid );
}
  
} // namespace mesh
