#include<stdexcept>
#include<iostream>

#include "Utility/Exceptions.hpp"
#include "EntityContainer.hpp"
#include "Entity.hpp"
#include "EntityStats.hpp"

namespace NSMesh { namespace NSEntity {

template<typename T>
EntityContainer<T>::EntityContainer( const NSMesh::Topology& p_topology ) : _topology( p_topology ) {
  // Create the entity stats class
  _entityStats = std::make_unique<EntityStats>();
}

    
template<typename T>
EntityContainer<T>::~EntityContainer() {
}


template<typename T>
void
EntityContainer<T>::AllocateSizeForEntity( size_t p_size ) {
  _entityVec.resize( p_size );
}


template<typename T>
void
EntityContainer<T>::SetEntityOwnershipIndexes( const NSMesh::Ownership& p_ownership, size_t p_size, size_t p_startIndex ) {
  NSExceptions::RuntimeError( _ownershipToEntityIndexesMap.find( p_ownership ) != _ownershipToEntityIndexesMap.end(), "You are setting the indexes for a given entity and ownership type that has been already set" );
  
  std::vector<std::pair<long unsigned, long unsigned> > tmp_vec;
  tmp_vec.push_back( std::make_pair( p_startIndex, p_startIndex + p_size ) );
  _ownershipToEntityIndexesMap.insert( std::make_pair( p_ownership, tmp_vec ) );
}


template<>    
void
EntityContainer<std::unique_ptr<NSMesh::NSEntity::Entity> >::AddEntityAtIndex( std::unique_ptr<NSMesh::NSEntity::Entity>& p_entity, unsigned p_gid, size_t p_index ) {
  NSExceptions::RuntimeError( p_index >= _entityVec.size(), "You are adding an entity to an entity container but you have not allocated enough memory to add that entity" );

  _entityVec[p_index] = std::move( p_entity );
  _gidToLidMap.insert( std::make_pair( p_gid, p_index ) );
}


template<>    
void
EntityContainer<NSMesh::NSEntity::Entity*>::AddEntityAtIndex( NSMesh::NSEntity::Entity*& p_entity, unsigned p_gid, size_t p_index ) {
  NSExceptions::RuntimeError( p_index >= _entityVec.size(), "You are adding an entity to an entity container but you have not allocated enough memory to add that entity" );
  
  _entityVec[p_index] = p_entity;
  _gidToLidMap.insert( std::make_pair( p_gid, p_index ) );
}

    
    
template<typename T>    
const std::vector<std::pair<long unsigned, long unsigned> >&
EntityContainer<T>::GetEntityIndexVec( const NSMesh::Ownership& p_ownership ) const {
  // Check if the entities you want to access are stored
  NSExceptions::RuntimeError( _ownershipToEntityIndexesMap.find( p_ownership ) == _ownershipToEntityIndexesMap.end(), "You are trying to access entities that you didn't store in the entities vector" );
  
  return _ownershipToEntityIndexesMap.find( p_ownership )->second;
}


template<typename T>    
const std::vector<T>&
EntityContainer<T>::GetEntityVec() const {
  return _entityVec;
}


template<typename T>    
const T&
EntityContainer<T>::GetEntity( const unsigned& p_gid ) const {
  // Get the lid
  const unsigned& lid = _gidToLidMap.find( p_gid )->second;
  
  return _entityVec[lid];
}


template<typename T>    
size_t
EntityContainer<T>::GetNumberOfEntity( const std::vector<NSMesh::Ownership>& p_ownership ) const {
  return _entityStats->GetNumberOfEntity( p_ownership );
}


template<typename T>    
void
EntityContainer<T>::SetNumberOfEntity( const NSMesh::Ownership& p_ownership, size_t p_number ) {
  _entityStats->SetNumberOfEntity( p_ownership, p_number );  
}


template<typename T>    
const size_t&
EntityContainer<T>::GetGlobalNumberOfEntity() const {
  return _entityStats->GetGlobalNumberOfEntity();
}


template<typename T>    
void
EntityContainer<T>::SetGlobalNumberOfEntity( const size_t& p_number ) {
  _entityStats->SetGlobalNumberOfEntity( p_number );
}
    

template<typename T>    
void
EntityContainer<T>::RebuildGIDToLIDMap() {
  // Clear the existing map
  _gidToLidMap.clear();

  // Loop over the entities vec
  for( size_t e=0; e<_entityVec.size(); ++e ) {
    // Get the entity
    const T& entity = _entityVec[e];

    // Add the gid to lid pair to the map
    _gidToLidMap.insert( std::make_pair( entity->GetGid(), entity->GetLid() ) );
  }
}

    
template<typename T>
bool
EntityContainer<T>::IsPresent( unsigned p_gid ) const {
  return ( _gidToLidMap.find( p_gid ) != _gidToLidMap.end() );
}
    
template class EntityContainer<std::unique_ptr<NSMesh::NSEntity::Entity> >;
template class EntityContainer<NSMesh::NSEntity::Entity*>;
    
} // namespace Entity
} // namespace Mesh
