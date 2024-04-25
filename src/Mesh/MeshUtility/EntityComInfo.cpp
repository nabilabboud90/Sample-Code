#include "Mesh/MeshUtility/EntityComInfo.hpp"
#include "Utility/Exceptions.hpp"

namespace NSMesh { namespace NSFields {

EntityComInfo::EntityComInfo() {
}


EntityComInfo::~EntityComInfo() {
}


EntityComInfo::EntityComInfo( const EntityComInfo& p_entityComInfo ) {
  _ownerProc = p_entityComInfo.GetOwningProcID();
  _sharingProcs = p_entityComInfo.GetSharingProcID();
  _universalProcs = p_entityComInfo.GetUniversalProcID();
}

    
EntityComInfo&
EntityComInfo::operator=( const EntityComInfo& p_entityComInfo ) {
  _ownerProc = p_entityComInfo.GetOwningProcID();
  _sharingProcs = p_entityComInfo.GetSharingProcID();
  _universalProcs = p_entityComInfo.GetUniversalProcID();

  return *this;
}

    
const int&
EntityComInfo::GetOwningProcID() const {
  return _ownerProc;
}


const std::vector<int>&
EntityComInfo::GetSharingProcID() const {
  return _sharingProcs;
}


const std::vector<int>&
EntityComInfo::GetUniversalProcID() const {
  return _universalProcs;
}


void
EntityComInfo::SetOwningProcID( int p_id ) {
  _ownerProc = p_id;
}


void
EntityComInfo::SetSharingProcID( const std::vector<int>& p_ids ) {
  _sharingProcs = p_ids;
}


void
EntityComInfo::SetUniversalProcID( const std::vector<int>& p_ids ) {
  _universalProcs = p_ids;
}


void
EntityComInfo::AddSharingProcID( int p_id ) {
  _sharingProcs.push_back( p_id );
}


void
EntityComInfo::AddUniversalProcID( int p_id ) {
  _universalProcs.push_back( p_id );
}
  
} // namespace MeshField
} // namespace Mesh
