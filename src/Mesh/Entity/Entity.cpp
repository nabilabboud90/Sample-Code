#include "Entity.hpp"
#include <iostream>

namespace NSMesh { namespace NSEntity {

Entity::Entity( const NSMesh::Topology& p_topology, unsigned p_gid, unsigned p_lid ) {
  _topology = p_topology;
  _gid = p_gid;
  _lid = p_lid;

  // Add the entity itself to the connectivity map
  // Class will return itself if the code asks about the connected entity
  _connectivity.resize( NSMesh::Topology::tp_LAST + 1, std::vector<const NSMesh::NSEntity::Entity*>() );
  _connectivity[ _topology ] = std::vector<const NSMesh::NSEntity::Entity*>( 1, this );
}


Entity::Entity( const Entity& p_e ) {
  _topology = p_e.GetTopology();
  _gid = p_e.GetGid();
  _lid = p_e.GetLid();
  
  // Add the entity itself to the connectivity map
  // Class will return itself if the code asks about the connected entity
  _connectivity.resize( NSMesh::Topology::tp_LAST + 1, std::vector<const NSMesh::NSEntity::Entity*>() );
  _connectivity[ _topology ] = std::vector<const NSMesh::NSEntity::Entity*>( 1, this );
}


void
Entity::SetGid( unsigned p_gid ) {
  _gid = p_gid;
}


void
Entity::SetLid( unsigned p_lid ) {
  _lid = p_lid;
}


void
Entity::SetTopology( const NSMesh::Topology& p_topology ) {
  _topology = p_topology;
}


unsigned
Entity::GetGid() const {
  return _gid;
}


const unsigned&
Entity::GetLid() const {
  return _lid;
}
    

NSMesh::Topology
Entity::GetTopology() const {
  return _topology;
}


void
Entity::AddConnectivity( const NSMesh::Topology& p_topology, const NSMesh::NSEntity::Entity* p_entity ) {
  if( _connectivity[ p_topology ].size() > 0 ) {
    std::vector<const NSMesh::NSEntity::Entity*>& entitiesVec = _connectivity[ p_topology ];

    // Check that the entity has not been added
    if( std::find( entitiesVec.begin(), entitiesVec.end(), p_entity ) == entitiesVec.end() ) {
      entitiesVec.push_back( p_entity );
    }
  }
  else {
    std::vector<const NSMesh::NSEntity::Entity*> entitiesVec( 1, p_entity );
    _connectivity[ p_topology ] = entitiesVec;
  }
}

    
int
Entity::GetNumConnectivity( const NSMesh::Topology& p_topology ) const {
  return _connectivity[ p_topology ].size();
}


const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >&
Entity::GetConnectivityInfo() const {
  return _connectivity;
}

} // namespace Entity
} // namespace Mesh
