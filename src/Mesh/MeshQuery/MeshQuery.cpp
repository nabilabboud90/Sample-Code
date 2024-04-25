#include "Mesh/MeshQuery/MeshQuery.hpp"
#include "Utility/Exceptions.hpp"
#include "Parser/InputParser/Datafile.hpp"

namespace NSMesh {

MeshQuery::MeshQuery() {
}


MeshQuery::~MeshQuery() {
}


void
MeshQuery::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  _name = p_sectionName;

  // Read the mesh query type
  std::string type = p_datafile.GetIfPresent<std::string>( p_sectionPath + "/" + p_sectionName + "/TYPE", "MESHQUERYR0" );

  if( type == "MESHQUERYR0" ) {
    _type = NSMesh::MeshQueryType::mqt_R0;
  }
  else if( type == "MESHQUERYR2" ) {
    _type = NSMesh::MeshQueryType::mqt_R2;
  }
  else {
    NSExceptions::InvalidArgument( true, "The mesh query type used is not defined in the code" );
  }
}


void
MeshQuery::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_mpiCom, int p_nDim ) {
  // Set the mpi communicator
  _mpiCom = p_mpiCom;
  
  // Set the dimension of the problem
  _ndim = p_nDim;
}


const std::map<NSMesh::Topology, std::vector<NSMesh::Topology> >&
MeshQuery::GetConnectivityInfoProtocol() const {
  return _connectivityInfoProtocol;
}


int
MeshQuery::GetSpatialDimension() const {
  return _ndim;
}


NSMesh::MeshQueryType
MeshQuery::GetType() const {
  return _type;
}

} // namespace MeshQuery
