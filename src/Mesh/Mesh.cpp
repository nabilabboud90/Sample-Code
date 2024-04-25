#include<iostream>

#include "Mesh.hpp"
#include "Entity/EntityContainer.hpp"
#include "Entity/Bucket.hpp"
#include "Mesh/MeshQuery/MeshQueryR2.hpp"
#include "Mesh/MeshQuery/MeshQueryR0.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "Communicator/MPI/MPIComWrapper.hpp"
#include "Utility/Exceptions.hpp"

namespace NSMesh {

Mesh::Mesh() {
}

  
void
Mesh::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  // Read the mesh file name
  _name = p_datafile.Get<std::string>( p_sectionPath + "/" + p_sectionName + "/NAME" );

  // Boolean specifying whether to create edges and faces on the whole mesh or just on parts of the mesh where boundary or interface conditions will be enforced
  _createAllSideEntities = p_datafile.GetIfPresent<bool>( p_sectionPath + "/" + p_sectionName + "/ALLSIDES", true );

  // Setup the mesh query class
  _setupMeshQuery( p_datafile, p_sectionPath, p_sectionName );
  
  // Read the mesh parts
  _setupMeshParts( p_datafile, p_sectionPath, p_sectionName );
}


void
Mesh::_setupMeshQuery( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  std::string meshQueryName = p_datafile.GetIfPresent<std::string>( p_sectionPath + "/" + p_sectionName + "/MESHQUERYNAME", "" );
  std::string meshQueryType = "MESHQUERYR0";
  if( meshQueryName != "" ) {
    meshQueryType = p_datafile.Get<std::string>( p_sectionPath + "/" + p_sectionName + "/" + meshQueryName + "/TYPE" );
  }
  _meshQuery.reset( MeshQuery::factory_mesh_query_Type::Create( meshQueryType ) );
  _meshQuery->Setup( p_datafile, p_sectionPath + "/" + p_sectionName, meshQueryName );
}

  
void
Mesh::_setupMeshParts( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  std::vector<std::string> mesh_parts_vec;
  p_datafile.Get( p_sectionPath + "/" + p_sectionName + "/LIST", mesh_parts_vec );
  
  // Check that all the different mesh parts have unique names
  std::vector<std::string>::iterator it_mesh_parts_name = std::unique( mesh_parts_vec.begin(), mesh_parts_vec.end() );
  bool was_mesh_parts_name_unique = ( it_mesh_parts_name == mesh_parts_vec.end() );
  NSExceptions::InvalidArgument( !was_mesh_parts_name_unique, "Mesh parts cannot have duplicate names" );

  // Loop over the different mesh parts and store the needed information
  for( size_t i=0; i<mesh_parts_vec.size(); ++i ) {
    // Read the type of the mesh part
    std::string type = p_datafile.Get<std::string>( p_sectionPath + "/" + p_sectionName + "/" + mesh_parts_vec[i] + "/TYPE" );
    if( type == "INTERIOR" ) {
      // Create the mesh part
      std::unique_ptr<NSMesh::InteriorMeshPart> interior_mesh_part_ptr;
      interior_mesh_part_ptr.reset( new NSMesh::InteriorMeshPart );
      
      // Setup the mesh part
      interior_mesh_part_ptr->Setup( p_datafile, p_sectionPath + "/" + p_sectionName, mesh_parts_vec[i] );
      
      // Add the mesh part to the mesh part name to mesh part map
      _nameToInteriorMeshPartMap.insert( std::make_pair( mesh_parts_vec[i], std::move( interior_mesh_part_ptr ) ) );
    }
    else if( type == "BOUNDARY" ) {
      // Create the mesh part
      std::unique_ptr<NSMesh::BoundaryMeshPart> boundary_mesh_part_ptr;
      boundary_mesh_part_ptr.reset( new NSMesh::BoundaryMeshPart );
      
      // Setup the mesh part
      boundary_mesh_part_ptr->Setup( p_datafile, p_sectionPath + "/" + p_sectionName, mesh_parts_vec[i] );
      
      // Add the mesh part to the mesh part name to mesh part map
      _nameToBoundaryMeshPartMap.insert( std::make_pair( mesh_parts_vec[i], std::move( boundary_mesh_part_ptr ) ) );
    }
    else {
      NSExceptions::InvalidArgument( true, "The type of the mesh part defined is not supported by the code" );
    }
  }
}


void
Mesh::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  // Set the mpi communicator
  _com = p_com;

  // Create the entity containers classes for the different entity types defined on the mesh (i.e. for vertices, elements, edges and faces)
  _initializeTopologyToEntityContainerMap();

  // Populate the entity container classes for each entity type defined on the mesh
  _populateTopologyToEntityContainerMap();

  // Build the connectivity relations
  // First initialize the mesh query class
  _meshQuery->Initialize( p_com, _ndim );

  // Use the connectivity information stored in the mesh query class to store the connectivities between the different mesh entities
  // This method below will store the connectivity information of the mesh in the Entity class
  // Note that the MeshQuery class specifies which connectivity information we want to store for each entity
  // That is the mesh query class specifies what connectivity information will be stored by the vertices, edges, faces, and elements
  _populateEntityConnectivityInformation();

  // Populate the mesh parts
  // Each MeshPart holds a vector of EntityContainer
  // Each index of this vector holds a pointer to an EntityContainer which holds information about the entities of a specific topology that are defined on the MeshPart
  _populateMeshPartsEntities();

  // Populate the communication information about the shared and universal entities
  _populateEntityComInfo();

  // Renumber the entities
  _renumberEntities();

  // Print a summary of the mesh data
  _printMeshSummary();
}


void
Mesh::_initializeTopologyToEntityContainerMap() {
  // Create the mesh vertices container
  auto vertex_container = std::make_unique<NSMesh::NSEntity::EntityContainer<std::unique_ptr<NSMesh::NSEntity::Entity> > >( NSMesh::Topology::tp_VERTEX );
  _topologyToEntityContainerMap.insert( std::make_pair( NSMesh::Topology::tp_VERTEX, std::move( vertex_container ) ) );

  // Create the mesh elements container
  auto element_container = std::make_unique<NSMesh::NSEntity::EntityContainer<std::unique_ptr<NSMesh::NSEntity::Entity> > >( NSMesh::Topology::tp_ELEMENT );
  _topologyToEntityContainerMap.insert( std::make_pair( NSMesh::Topology::tp_ELEMENT, std::move( element_container ) ) );

  // Create the mesh edges container
  if( _ndim > 1 ) {
    auto edge_container = std::make_unique<NSMesh::NSEntity::EntityContainer<std::unique_ptr<NSMesh::NSEntity::Entity> > >( NSMesh::Topology::tp_EDGE );    
    _topologyToEntityContainerMap.insert( std::make_pair( NSMesh::Topology::tp_EDGE, std::move( edge_container ) ) );
  }

  // Create the mesh faces container
  if( _ndim == 3 ) {
    auto face_container = std::make_unique<NSMesh::NSEntity::EntityContainer<std::unique_ptr<NSMesh::NSEntity::Entity> > >( NSMesh::Topology::tp_FACE );
    _topologyToEntityContainerMap.insert( std::make_pair( NSMesh::Topology::tp_FACE, std::move( face_container ) ) );
  }
}


void FillBucket( const NSMesh::Mesh* const p_mesh, NSEntity::Bucket<const NSEntity::Entity*>& p_bucket, const NSMesh::Topology& p_topology,
		 const std::vector<unsigned>& p_entityGids, unsigned& p_entityCnt  ) {
  // Loop over the capacity of the bucket
  std::vector<const NSMesh::NSEntity::Entity*> entities( p_bucket.GetSize(), nullptr );
  for( size_t e=0; e<p_bucket.GetSize(); ++e ) {
    // Get the entity
    const auto& entity = p_mesh->GetEntity( p_topology, p_entityGids[ p_entityCnt ] );
    
    // Set the pointer to the entity in the bucket
    entities[ e ] = &entity;
    
      // Increment the entity counter
    p_entityCnt++;
  }
  
    // Set the vector holding the pointers to the entities to be referenced by the bucket
  p_bucket.SetValues( entities );
}

  
void
Mesh::GetBucketsByGid( std::vector<NSEntity::Bucket<const NSEntity::Entity*> >& p_buckets, const std::vector<unsigned>& p_entityGids,
		       const NSMesh::Topology& p_topology, unsigned p_capacity ) const {
  // Check that the buckets that are being passed are empty
  NSExceptions::LogicError( p_buckets.size() > 0, "You are trying to fill a vector of bucket that already have some bucket in it" );

  // Get the number of buckets
  unsigned remainder = p_entityGids.size() % p_capacity;
  unsigned total_num_buckets = ( p_entityGids.size() - remainder ) / p_capacity + 1;
  unsigned num_buckets_with_full_capacity = total_num_buckets - 1;
  
  // Resize the p_buckets vector
  p_buckets.resize( total_num_buckets );
  
  // Fill the buckets with full capacity
  unsigned entity_cnt = 0;
  for( size_t b=0; b<num_buckets_with_full_capacity; ++b ) {
    // Get a reference to the bucket
    auto& bucket = p_buckets[b];

    // Set the size of the bucket
    bucket.Resize( p_capacity );

    // Fill the bucket
    FillBucket( this, bucket, p_topology, p_entityGids, entity_cnt );
  }

  // Fill any buckets whose size is less than the full capacity
  if( remainder != 0 ) {
    // Get a reference to the bucket
    auto& bucket = p_buckets[ num_buckets_with_full_capacity ];
    
    // Set the size
    bucket.Resize( remainder );

    // Fill the bucket
    FillBucket( this, bucket, p_topology, p_entityGids, entity_cnt );
  }
}

  
void
Mesh::GetBuckets( std::vector<NSMesh::NSEntity::Bucket<std::unique_ptr<NSMesh::NSEntity::Entity> > >& p_buckets,
		  const std::vector<NSMesh::Ownership>& p_ownership, const NSMesh::Topology& p_topology, unsigned p_capacity ) const {
  // First get the entity container for specific entity topology (i.e, vertex, element, face, or edge)
  NSExceptions::RuntimeError( _topologyToEntityContainerMap.find( p_topology ) == _topologyToEntityContainerMap.end(), "You are trying to gather entities of a certain topology that you have not stored" );

  // Second check that the buckets that are being passed are empty
  NSExceptions::LogicError( p_buckets.size() > 0, "You are trying to fill a vector of bucket that already have some bucket in it" );

  const auto& entity_container = _topologyToEntityContainerMap.find( p_topology )->second;

  // Third get the vector specifying the indexes of the entites of the specified topology and ownership (i.e, owned vertices, shared elements, etc, or a union of owned and shared elements for example) wihtin the entity vector
  std::vector<std::pair<long unsigned, long unsigned> > indexes_vec;
  for( size_t i=0; i<p_ownership.size(); ++i ) {
    // Get the index vec for a given ownership
    const std::vector<std::pair<long unsigned, long unsigned> >& single_ownership_index_vec = entity_container->GetEntityIndexVec( p_ownership[i] );

    for( size_t o=0; o<single_ownership_index_vec.size(); ++o ) {
      // Push back the indexes for a given ownership to the indexes vec
      indexes_vec.push_back( single_ownership_index_vec[o] );
    }
  }

  // Get the entity vector
  const std::vector<std::unique_ptr<NSMesh::NSEntity::Entity> >& entity_vec = entity_container->GetEntityVec();

  // Populate the buckets vector
  _populateBuckets<std::unique_ptr<NSMesh::NSEntity::Entity>, std::unique_ptr<NSMesh::NSEntity::Entity> >( p_buckets, indexes_vec, entity_vec, p_capacity );
}


void
Mesh::GetBucketsByBoundaryPart( std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> >& p_buckets, const std::vector<NSMesh::Ownership>& p_ownership,
				const NSMesh::Topology& p_topology, const std::string& p_boundaryPartName, const std::string& p_interiorPartName, unsigned p_capacity ) const {
  // Check that the topology is either element or edge/face if the problem is 2D/3D
  NSExceptions::RuntimeError( p_topology != NSMesh::Topology::tp_ELEMENT && p_topology != GetClosureTopology(), "Boundary mesh parts only store information about Elements and Closure entities for now" );
  
  // Get the entity container for specific entity topology (i.e, vertex, element, face, or edge)
  NSExceptions::RuntimeError( _topologyToEntityContainerMap.find( p_topology ) == _topologyToEntityContainerMap.end(), "You are trying to gather entities of a certain topology that you have not stored" );
  
  // Check that the buckets that are being passed are empty
  NSExceptions::LogicError( p_buckets.size() > 0, "You are trying to fill a vector of bucket that already have some bucket in it" );
  
  // Check that the mesh part is defind
  const auto it_mesh_part = _nameToBoundaryMeshPartMap.find( p_boundaryPartName );
  NSExceptions::LogicError( it_mesh_part == _nameToBoundaryMeshPartMap.end(), "You are trying to get the entities from a boundary mesh part that is not defined" );

  // Get the mesh part
  const std::unique_ptr<NSMesh::BoundaryMeshPart>& mesh_part = it_mesh_part->second;
  
  // Second get the vector specifying the indexes of the entites of the specified topology and ownership (i.e, owned vertices, shared elements, etc, or a union of owned and shared elements for example) wihtin the entity vector
  std::vector<std::pair<long unsigned, long unsigned> > indexes_vec;
  for( size_t i=0; i<p_ownership.size(); ++i ) {
    // Get the index vec for a given ownership
    const auto& single_ownership_index_vec = mesh_part->GetEntityIndexVec( p_topology, p_interiorPartName, p_ownership[i] );

    for( size_t o=0; o<single_ownership_index_vec.size(); ++o ) {
      // Push back the indexes for a given ownership to the indexes vec
      indexes_vec.push_back( single_ownership_index_vec[o] );
    }
  }

  // Get the entity vector
  const auto& entity_vec = mesh_part->GetEntityVec( p_topology, p_interiorPartName );
  
  // Populate the buckets vector
  _populateBuckets<const NSMesh::NSEntity::Entity*, NSMesh::NSEntity::Entity*>( p_buckets, indexes_vec, entity_vec, p_capacity );
}
  

void
Mesh::GetBucketsByInteriorPart( std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> >& p_buckets, const std::vector<NSMesh::Ownership>& p_ownership, const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, unsigned p_capacity ) const {
  // Get the entity container for specific entity topology (i.e, vertex, element, face, or edge)
  NSExceptions::RuntimeError( _topologyToEntityContainerMap.find( p_topology ) == _topologyToEntityContainerMap.end(), "You are trying to gather entities of a certain topology that you have not stored" );
  
  // Check that the buckets that are being passed are empty
  NSExceptions::LogicError( p_buckets.size() > 0, "You are trying to fill a vector of bucket that already have some bucket in it" );
  
  // Check that the mesh part is defind
  const auto it_mesh_part = _nameToInteriorMeshPartMap.find( p_interiorPartName );
  NSExceptions::LogicError( it_mesh_part == _nameToInteriorMeshPartMap.end(), "You are trying to get the entities from an interior mesh part that is not defined" );

  // Get the mesh part
  const std::unique_ptr<NSMesh::InteriorMeshPart>& mesh_part = it_mesh_part->second;
  
  // Second get the vector specifying the indexes of the entites of the specified topology and ownership (i.e, owned vertices, shared elements, etc, or a union of owned and shared elements for example) wihtin the entity vector
  std::vector<std::pair<long unsigned, long unsigned> > indexes_vec;
  for( size_t i=0; i<p_ownership.size(); ++i ) {
    // Get the index vec for a given ownership
    const auto& single_ownership_index_vec = mesh_part->GetEntityIndexVec( p_topology, p_ownership[i] );

    for( size_t o=0; o<single_ownership_index_vec.size(); ++o ) {
      // Push back the indexes for a given ownership to the indexes vec
      indexes_vec.push_back( single_ownership_index_vec[o] );
    }
  }

  // Get the entity vector
  const auto& entity_vec = mesh_part->GetEntityVec( p_topology );
  
  // Populate the buckets vector
  _populateBuckets<const NSMesh::NSEntity::Entity*, NSMesh::NSEntity::Entity*>( p_buckets, indexes_vec, entity_vec, p_capacity );
}


template<typename T1, typename T2>
void
Mesh::_populateBuckets( std::vector<NSMesh::NSEntity::Bucket<T1> >& p_buckets,
			const std::vector<std::pair<long unsigned, long unsigned> >& p_indexesVec,
			const std::vector<T2>& p_entityVec, unsigned p_capacity ) const {
  // Loop a first time over the vector of pair of indexes to figure out the total number of buckets
  unsigned numBuckets = 0;
  for( size_t i=0; i<p_indexesVec.size(); ++i ) {
    // Get the size of the vector of buckets
    const std::pair<long unsigned, long unsigned>& index_pair = p_indexesVec[i];
    unsigned size = std::get<1>( index_pair ) - std::get<0>( index_pair );
    unsigned remainder = size%p_capacity;
    numBuckets += ( size - remainder ) / p_capacity;
    
    // Add a trailing bucket if any
    numBuckets += ( (remainder>0) ? 1 : 0 );
  }
  
  // Resize the buckets vector
  p_buckets.resize( numBuckets );
  
  // Loop a second time to set the information for each bucket in the buckets vector
  unsigned bucket_cnt = 0;
  numBuckets = 0;
  for( size_t i=0; i<p_indexesVec.size(); ++i ) {
    // Get the size of the vector of buckets
    const std::pair<long unsigned, long unsigned>& index_pair = p_indexesVec[i];

    // Set the starting index
    unsigned starting_index = std::get<0>( index_pair );

    // Get the number of buckets
    unsigned size = std::get<1>( index_pair ) - std::get<0>( index_pair );
    unsigned remainder = size%p_capacity;    
    numBuckets += ( size - remainder ) / p_capacity;

    // Loop over each bucket and set its starting index in the entites vec
    // Note that numBuckets represents all the full buckets excluding any bucket that will hold the remainder of the total number of entities
    for( size_t j=bucket_cnt; j<numBuckets; ++j ) {
      p_buckets[j].Resize( p_capacity );
      auto& data_ptr = ( p_entityVec[ starting_index ] );
      p_buckets[j].SetStartingIndex( &data_ptr );
      starting_index += p_capacity;
    }
    
    // Check if there is any trailing bucket which has a size less than capacity and specify its starting index in the entities vec
    if( remainder > 0 ) {
      p_buckets[numBuckets].Resize( remainder );
      auto& data_ptr = ( p_entityVec[ starting_index ] );
      p_buckets[numBuckets].SetStartingIndex( &data_ptr );
      numBuckets++;
    }
    bucket_cnt = numBuckets;
  }
}
  

const NSMesh::NSEntity::Entity&
Mesh::GetEntity( const NSMesh::Topology& p_topology, const unsigned& p_gid ) const {
  // Get the entity container corresponding to the given entity topology we are asking for
  const auto& entity_container_ptr = _topologyToEntityContainerMap.find( p_topology )->second;
  return *entity_container_ptr->GetEntity( p_gid );
}


NSMesh::Topology
Mesh::GetEntityTopology( const unsigned& p_gid ) const {  
  for( const auto& info : _topologyToEntityContainerMap ) {
    // Get the pointer to the entity container
    const auto& entity_container = info.second;

    if( entity_container->IsPresent( p_gid ) ) {
      return info.first;
    }
  }

  NSExceptions::InvalidArgument( true, "There is no entity in the mesh with gid " + p_gid );
  return NSMesh::Topology::tp_VERTEX;
}
  

void
Mesh::Begin( const NSMesh::Topology& p_toTopology, const NSMesh::NSEntity::Entity& p_entity, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVector ) const {
  // Ask the entity to return its connectivity map
  const auto& connectivityInfo = p_entity.GetConnectivityInfo();

  // Pass the connectivity map to the mesh query class which has the particular functions that recover the connected entities of the given topology
  _meshQuery->Begin( p_entity.GetTopology(), p_toTopology, connectivityInfo, p_connectedEntitiesVector );
}


const std::string&
Mesh::GetName() const {
  return _name;
}


unsigned
Mesh::GetNumberOfEntity( const NSMesh::Topology& p_topology, const std::vector<NSMesh::Ownership>& p_ownership ) const {
  NSExceptions::LogicError( _topologyToEntityContainerMap.find( p_topology ) == _topologyToEntityContainerMap.end(), "You are asking for the number of entities of an entity which is not on this mesh" );

  return _topologyToEntityContainerMap.find( p_topology )->second->GetNumberOfEntity( p_ownership );
}


unsigned
Mesh::GetNumberOfEntityPerParts( const NSMesh::Topology& p_topology, const std::vector<NSMesh::Ownership>& p_ownership, const std::vector<std::string>& p_parts ) const {      
  unsigned total_num_entities = 0;
  for( size_t p=0; p<p_parts.size(); ++p ) {
    NSExceptions::LogicError( _nameToInteriorMeshPartMap.find( p_parts[p] ) == _nameToInteriorMeshPartMap.end(), "The GetNumberOfEntityPerParts is currently only supported for interior mesh parts" );  
    NSExceptions::LogicError( _topologyToEntityContainerMap.find( p_topology ) == _topologyToEntityContainerMap.end(), "You are asking for the number of entities of an entity which is not on this mesh" );

    total_num_entities += _nameToInteriorMeshPartMap.find( p_parts[p] )->second->GetNumberOfEntity( p_topology, p_ownership );
  }
  
  return total_num_entities;
}
  

const size_t&
Mesh::GetGlobalNumberOfEntity( const NSMesh::Topology& p_topology ) const {
  NSExceptions::LogicError( _topologyToEntityContainerMap.find( p_topology ) == _topologyToEntityContainerMap.end(), "You are asking for the global number of entities of an entity which is not on this mesh" );

  return _topologyToEntityContainerMap.find( p_topology )->second->GetGlobalNumberOfEntity();
}


int
Mesh::SpatialDimension() const {
  return _ndim;
}


const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& 
Mesh::GetMPIComWrapper() const {
  return _com;
}


const std::map<NSMesh::Topology, std::map<unsigned, NSMesh::NSFields::EntityComInfo> >&
Mesh::GetProcInterfaceEntityInfo() const {
  return _topoToEntityGidProcMap;
}


void
Mesh::_printMeshSummary() const {
  int rank = _com->GetRank();
  if( rank == 0 ) {
    // Get the element stats
    size_t totalNumElements = _topologyToEntityContainerMap.find( NSMesh::Topology::tp_ELEMENT )->second->GetGlobalNumberOfEntity();

    // Get the vertices stats
    size_t totalNumVertices = _topologyToEntityContainerMap.find( NSMesh::Topology::tp_VERTEX )->second->GetGlobalNumberOfEntity();

    // Get the edges stats
    size_t totalNumEdges = _topologyToEntityContainerMap.find( NSMesh::Topology::tp_EDGE )->second->GetGlobalNumberOfEntity();

    // Get the faces stats
    size_t totalNumFaces = 0;
    if( _ndim == 3 ) {
      totalNumFaces = _topologyToEntityContainerMap.find( NSMesh::Topology::tp_FACE )->second->GetGlobalNumberOfEntity();
    }

    std::cout << "************************************************************" << std::endl;
    std::cout << "*******************  PRINTING MESH INFO ********************" << std::endl;
    std::cout << "************************************************************" << std::endl;
    std::cout << "Total number of vertices: " << totalNumVertices << std::endl;
    std::cout << "Total number of edges: " << totalNumEdges << std::endl;
    std::cout << "Total number of faces: " << totalNumFaces << std::endl;
    std::cout << "Total number of elements: " << totalNumElements << std::endl;
    std::cout << "************************************************************" << std::endl;
  }
}


bool
Mesh::IsMeshPartDefined( const std::string& p_name ) const {
  if( _nameToInteriorMeshPartMap.find( p_name ) != _nameToInteriorMeshPartMap.end() ) {
    return true;
  }
  if( _nameToBoundaryMeshPartMap.find( p_name ) != _nameToBoundaryMeshPartMap.end() ) {
    return true;
  }

  return false;
}


NSMesh::Topology
Mesh::GetClosureTopology() const {
  if( _ndim == 2 ) {
    return NSMesh::Topology::tp_EDGE;
  }
  else if( _ndim == 3 ) {
    return NSMesh::Topology::tp_FACE;
  }
  else {
    NSExceptions::LogicError( true, "The closure topology is not defined for the spatial dimension of the problem being simulated" );
    return NSMesh::Topology::tp_VERTEX;
  }
}


bool
Mesh::IsOnProcessorInterface( unsigned p_gid ) const {
  // Loop over the data structure holding information about the interface entitites
  for( const auto& entity_interface_info : _topoToEntityGidProcMap ) {
    const auto& proc_entities_with_specific_topology = entity_interface_info.second;
  
    // Check if the gid of the entity is stored in the above data structure
    if( proc_entities_with_specific_topology.find( p_gid ) != proc_entities_with_specific_topology.end() ) {
      return true;
    }
  }

  return false;
}


bool
Mesh::IsOwned( unsigned p_gid ) const {
  // Loop over the data structure holding information about the interface entitites
  for( const auto& entity_interface_info : _topoToEntityGidProcMap ) {
    const auto& proc_entities_with_specific_topology = entity_interface_info.second;
    
    // Check if the gid of the entity is stored in the above data structure
    auto it = proc_entities_with_specific_topology.find( p_gid );
    if( it != proc_entities_with_specific_topology.end() ) {
      // Get the rank of the processor
      auto rank = _com->GetRank();

      // Check if the rank is equal to the rank of the owning processor
      if( it->second.GetOwningProcID() == rank ) {
	return true;
      }
      else {
	return false;
      }
    }
  }
  
  return true;  
}
  

unsigned
Mesh::GetOwningProcessor( unsigned p_gid ) const {
  // Loop over the data structure holding information about the interface entitites
  for( const auto& entity_interface_info : _topoToEntityGidProcMap ) {
    const auto& proc_entities_with_specific_topology = entity_interface_info.second;
    
    // Check if the gid of the entity is stored in the above data structure
    auto it = proc_entities_with_specific_topology.find( p_gid );
    if( it != proc_entities_with_specific_topology.end() ) {
      return it->second.GetOwningProcID();
    }
  }

  return _com->GetRank();
}
  
  
NSMesh::NSEntity::Entity&
Mesh::_getEntity( const NSMesh::Topology& p_topology, const unsigned& p_gid ) {
  // Get the entity container corresponding to the given entity topology we are asking for
  const auto& entity_container_ptr = _topologyToEntityContainerMap.find( p_topology )->second;  
  return *entity_container_ptr->GetEntity( p_gid );
}


const std::map<unsigned, NSMesh::NSFields::EntityComInfo>&
Mesh::GetProcessorInterfaceEntities( const NSMesh::Topology& p_topology ) const {
  NSExceptions::RuntimeError( _topoToEntityGidProcMap.find( p_topology ) == _topoToEntityGidProcMap.end(), "Trying to get information fomr _topoToEntityGidProcMap for a topology that was registered in the map" );
  return _topoToEntityGidProcMap.find( p_topology )->second;
}


std::vector<std::string>
Mesh::GetAllMeshPartsName( const NSMesh::Topology& p_topology ) const {
  if( p_topology == NSMesh::Topology::tp_ELEMENT ) {
    std::vector<std::string> result;
    for( const auto& part : _nameToInteriorMeshPartMap ) {
      // Get the mesh part name
      std::string part_name = part.first;
      result.push_back( part_name );
    }

    return result;
  }
  if( p_topology == GetClosureTopology() ) {
    std::vector<std::string> result;
    for( const auto& part : _nameToBoundaryMeshPartMap ) {
      // Get the mesh part name
      std::string part_name = part.first;
      result.push_back( part_name );
    }

    return result;
  }

  return std::vector<std::string>();
}  


NSMesh::ElementTopology
Mesh::GetPartElementTopology( const std::string& p_partName ) const {
  if( _nameToInteriorMeshPartMap.find( p_partName ) != _nameToInteriorMeshPartMap.end() ) {
    return _nameToInteriorMeshPartMap.find( p_partName )->second->GetMeshElementTopology();    
  }
  else if( _nameToBoundaryMeshPartMap.find( p_partName ) != _nameToBoundaryMeshPartMap.end() ) {
    return _nameToBoundaryMeshPartMap.find( p_partName )->second->GetMeshElementTopology();    
  }
  else {
    NSExceptions::InvalidArgument( true, "The part whose element topology is asked for is not defined" );
    return NSMesh::ElementTopology::et_UNDEFINED;
  }
}


int
Mesh::GetNumberOfVerticesPerElement( const std::string& p_partName ) const {
  if( _nameToInteriorMeshPartMap.find( p_partName ) != _nameToInteriorMeshPartMap.end() ) {
    return _nameToInteriorMeshPartMap.find( p_partName )->second->GetNumberOfVerticesPerElement();
  }
  else if( _nameToBoundaryMeshPartMap.find( p_partName ) != _nameToBoundaryMeshPartMap.end() ) {
    return _nameToBoundaryMeshPartMap.find( p_partName )->second->GetNumberOfVerticesPerElement();
  }
  else {
    NSExceptions::InvalidArgument( true, "The part whose number of mesh vertices per element is asked for is not defined" );
    return 0;
  }
}


const NSMesh::BoundaryMeshPart*
Mesh::GetBoundaryMeshPart( const std::string& p_partName ) const {
  NSExceptions::RuntimeError( _nameToBoundaryMeshPartMap.find( p_partName ) == _nameToBoundaryMeshPartMap.end(), "The boundary mesh part asked for in GetBoundaryMeshPart is not defined" );
  return _nameToBoundaryMeshPartMap.find( p_partName )->second.get();
}


const NSMesh::InteriorMeshPart*
Mesh::GetInteriorMeshPart( const std::string& p_partName ) const {
  NSExceptions::RuntimeError( _nameToInteriorMeshPartMap.find( p_partName ) == _nameToInteriorMeshPartMap.end(), "The interior mesh part asked for in GetInteriorMeshPart is not defined" );
  return _nameToInteriorMeshPartMap.find( p_partName )->second.get();
}


void
Mesh::SortByTopology( const std::vector<unsigned>& p_entitiesGid, std::vector<std::vector<unsigned> >& p_sortedEntitiesGid ) const {
  // Resize the sortedEntitiesGid vector to the size of all the entity topologies supported in the code
  p_sortedEntitiesGid.resize( 4, std::vector<unsigned>() );

  // Reserve some memory for each vector of a specific entity topology
  for( auto& entity_gid : p_sortedEntitiesGid ) {
    entity_gid.reserve( p_entitiesGid.size() );
  }
    
  // Loop over the entities gid vector
  for( const auto& gid : p_entitiesGid ) {
    // Get the topology of the entity
    auto topology = GetEntityTopology( gid );

    // Store the entity gid in the corresponding location
    p_sortedEntitiesGid[ topology ].push_back( gid );
  }

  // Remove any non used allocated memory
  for( auto& entity_gid : p_sortedEntitiesGid ) {
    entity_gid.shrink_to_fit();
  }
}
  
} // namespace Mesh
