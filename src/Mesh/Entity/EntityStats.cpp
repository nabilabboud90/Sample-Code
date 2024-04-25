#include<iostream>

#include "Utility/Exceptions.hpp"
#include "Mesh/Entity/EntityStats.hpp"

namespace NSMesh { namespace NSEntity {

EntityStats::EntityStats() {
}


EntityStats::~EntityStats() {
}


size_t
EntityStats::GetNumberOfEntity( const std::vector<NSMesh::Ownership>& p_ownership ) const {
  size_t result = 0.0;
  for( size_t i=0; i<p_ownership.size(); ++i ) {
    NSExceptions::LogicError( _ownershipToEntityNumberMap.find( p_ownership[i] ) == _ownershipToEntityNumberMap.end(), "Trying to get the number of entities of a given topology and ownership that was not stored" );

    result += _ownershipToEntityNumberMap.find( p_ownership[i] )->second;
  }

  return result;
}


void
EntityStats::SetNumberOfEntity( const NSMesh::Ownership& p_ownership, size_t p_number ) {
  _ownershipToEntityNumberMap.insert( std::make_pair( p_ownership, p_number ) );
}

    
const size_t&
EntityStats::GetGlobalNumberOfEntity() const {
  return _globalNumberEntity;
}


void
EntityStats::SetGlobalNumberOfEntity( const size_t& p_number ) {
  _globalNumberEntity = p_number;
}

} // namespace Entity
} // namespace Mesh
