#include<iostream>

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <Shards_CellTopology.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/FillMesh.hpp>
#pragma GCC diagnostic pop

#include "STKMesh.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "Communicator/MPI/MPIComWrapper.hpp"
#include "Mesh/Entity/EntityContainer.hpp"
#include "Mesh/Entity/Bucket.hpp"
#include "Mesh/Entity/Entity.hpp"
#include "Mesh/MeshQuery/MeshQueryR2.hpp"
#include "Mesh/MeshPart/InteriorMeshPart.hpp"
#include "Mesh/MeshPart/BoundaryMeshPart.hpp"
#include "Mesh/Entity/EntityStats.hpp"

namespace NSMesh { namespace NSSTKMesh {

Mesh* createSTKMesh() {
  return new STKMesh;
}


STKMesh::STKMesh() {
  _topologyToSTKMeshTopology.insert( std::make_pair( NSMesh::Topology::tp_VERTEX, stk::topology::NODE_RANK ) );
  _topologyToSTKMeshTopology.insert( std::make_pair( NSMesh::Topology::tp_EDGE, stk::topology::EDGE_RANK ) );
  _topologyToSTKMeshTopology.insert( std::make_pair( NSMesh::Topology::tp_FACE, stk::topology::FACE_RANK ) );
  _topologyToSTKMeshTopology.insert( std::make_pair( NSMesh::Topology::tp_ELEMENT, stk::topology::ELEMENT_RANK ) );
}


void
STKMesh::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  Mesh::Setup( p_datafile, p_sectionPath, p_sectionName );
}


void
STKMesh::_initializeMeshParts() {
  // Initialize the boundary mesh parts
  for( auto& boundary_part : _nameToBoundaryMeshPartMap ) {
    boundary_part.second->Initialize( _ndim );
  }
}


void
STKMesh::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  // Use STK IO to populate a STK Mesh
  _meshReader = std::make_unique<stk::io::StkMeshIoBroker>( p_com->GetMPICom() );

  stk::mesh::MeshBuilder builder(p_com->GetMPICom());
  builder.set_aura_option(stk::mesh::BulkData::AUTO_AURA);
  _bulkData = builder.create();

  // Populate bulk data
  _meshReader->set_bulk_data(*_bulkData);
  _meshReader->property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
  _meshReader->add_mesh_database(_name, stk::io::READ_MESH);
  _meshReader->create_input_mesh();
  _meshReader->populate_bulk_data();
  
  // Populate the meta data from STK IO
  _metaData = &_meshReader->meta_data();
  
  // Get the dimension of the problem
  _ndim = _metaData->spatial_dimension();

  // Initialize the mesh parts
  _initializeMeshParts();
  
  // Create the stk mesh parts equivalent to the mesh parts stored in 
  _createSTKMeshParts( p_com );

  // Create the faces and edges in the mesh
  if( _ndim > 1 ) {
    _createSTKEdges();
  }

  if( _ndim == 3 ) {
    _createSTKFaces();
  }

  // Call the parent's initialize function
  Mesh::Initialize( p_com );
}


void
STKMesh::_createSTKEdges() {
  if( _createAllSideEntities ) {
    stk::mesh::create_edges( *_bulkData );
  }
  else {
    // Loop over the stk parts
    for( const auto& it : _partNameToStkMeshPart ) {
      if( it.second->topology() == stk::topology::EDGE_RANK || it.second->topology() == stk::topology::FACE_RANK ) {
	stk::mesh::create_edges( *_bulkData, *(it.second) );
      }
    }
  }
}


void
STKMesh::_createSTKFaces() {
  if( _createAllSideEntities ) {
    stk::mesh::create_faces( *_bulkData, true );
  }
  else {
    // Loop over the stk parts
    for( const auto& it : _partNameToStkMeshPart ) {
      if( it.second->topology() == stk::topology::FACE_RANK ) {
	stk::mesh::create_faces( *_bulkData, *(it.second), true );
      }
    }
  }
}


void
STKMesh::_renumberEntities() {
  // Renumber entities
  // The renumbering ensures that each entity (vertex, element, face and edge) has a unique gid across processors
  // This means that if a vertex has a global id of 1 then no other entity of any topology can have this global id across all the processors
  size_t offset = 0;
  if( _renumberMeshEntities ) {
    std::map<unsigned, unsigned> oldToNewVerticesGIDMap;
    _renumberEntitiesImpl( NSMesh::Topology::tp_VERTEX, stk::topology::NODE_RANK, oldToNewVerticesGIDMap, offset );
    _communicateNewGIDs( NSMesh::Topology::tp_VERTEX, oldToNewVerticesGIDMap );
    offset += GetGlobalNumberOfEntity( NSMesh::Topology::tp_VERTEX );
    
    std::map<unsigned, unsigned> oldToNewElementsGIDMap;
    _renumberEntitiesImpl( NSMesh::Topology::tp_ELEMENT, stk::topology::ELEMENT_RANK, oldToNewElementsGIDMap, offset );
    _communicateNewGIDs( NSMesh::Topology::tp_ELEMENT, oldToNewElementsGIDMap );
    offset += GetGlobalNumberOfEntity( NSMesh::Topology::tp_ELEMENT );
    
    if( _ndim > 1 ) {
      std::map<unsigned, unsigned> oldToNewEdgesGIDMap;
      _renumberEntitiesImpl( NSMesh::Topology::tp_EDGE, stk::topology::EDGE_RANK, oldToNewEdgesGIDMap, offset );
      _communicateNewGIDs( NSMesh::Topology::tp_EDGE, oldToNewEdgesGIDMap );
      offset += GetGlobalNumberOfEntity( NSMesh::Topology::tp_EDGE );
    }

    if( _ndim == 3 ) {
      std::map<unsigned, unsigned> oldToNewFacesGIDMap;
      _renumberEntitiesImpl( NSMesh::Topology::tp_FACE, stk::topology::FACE_RANK, oldToNewFacesGIDMap, offset );
      _communicateNewGIDs( NSMesh::Topology::tp_FACE, oldToNewFacesGIDMap );
      offset += GetGlobalNumberOfEntity( NSMesh::Topology::tp_FACE );
    }
  }
}


void
STKMesh::_populateTopologyToEntityContainerMap() {
  // Store the information about the vertices of the mesh in the EntityContainer
  // This includes the vertex gids as well as the elements, edges, and faces to which they are connected
  _populateEntityInfo( stk::topology::NODE_RANK, NSMesh::Topology::tp_VERTEX );

  // Store the information about the elements of the mesh in the EntityContainer
  _populateEntityInfo( stk::topology::ELEMENT_RANK, NSMesh::Topology::tp_ELEMENT );;

  // Store the information about the edges of the mesh in the EntityContainer
  if( _ndim > 1 ) {
    _populateEntityInfo( stk::topology::EDGE_RANK, NSMesh::Topology::tp_EDGE );;
  }

  // Store the information about the faces of the mesh in the EntityContainer
  if( _ndim == 3 ) {
    _populateEntityInfo( stk::topology::FACE_RANK, NSMesh::Topology::tp_FACE );;
  }
}


void
STKMesh::_populateEntityInfo( const stk::mesh::EntityRank& p_rank, const NSMesh::Topology& p_topology ) {
  // Get the total number of entities that this processor has access to ( owned + shared + universal )
  size_t total_num_owned_entities = _getTotalNumEntities( p_rank, _metaData->locally_owned_part() );
  size_t total_num_shared_entities = _getTotalNumEntities( p_rank, _metaData->globally_shared_part() & !_metaData->locally_owned_part() );
  size_t total_num_universal_entities = _getTotalNumEntities( p_rank, _metaData->universal_part() &(!_metaData->locally_owned_part()) & (!_metaData->globally_shared_part()) );
  size_t total_num_entities = total_num_owned_entities + total_num_shared_entities + total_num_universal_entities;

  // Set the entity stats
  _topologyToEntityContainerMap.find( p_topology )->second->SetNumberOfEntity( NSMesh::Ownership::o_OWNED, total_num_owned_entities );
  _topologyToEntityContainerMap.find( p_topology )->second->SetNumberOfEntity( NSMesh::Ownership::o_SHARED, total_num_shared_entities );
  _topologyToEntityContainerMap.find( p_topology )->second->SetNumberOfEntity( NSMesh::Ownership::o_UNIVERSAL, total_num_universal_entities );
  
  // Get and set the total number of this entity on the whole mesh
  size_t global_num_entities = 0;
  _com->AllReduce( &total_num_owned_entities, &global_num_entities, 1, my_MPI_SIZE_T, MPI_SUM );
  _topologyToEntityContainerMap.find( p_topology )->second->SetGlobalNumberOfEntity( global_num_entities );

  // Allocate memory for the entities vector which will hold all the entities that his process has access to
  _topologyToEntityContainerMap.find( p_topology )->second->AllocateSizeForEntity( total_num_entities );
  
  // Set the starting and ending index of the owned entities within the vector that holds all the entities (owned, shared, and universal)
  _topologyToEntityContainerMap.find( p_topology )->second->SetEntityOwnershipIndexes( NSMesh::Ownership::o_OWNED, total_num_owned_entities, 0 );

  // Set the starting and ending index of the shared entities within the vector that holds all the entities (owned, shared, and universal)
  _topologyToEntityContainerMap.find( p_topology )->second->SetEntityOwnershipIndexes( NSMesh::Ownership::o_SHARED, total_num_shared_entities, total_num_owned_entities );

  // Set the starting and ending index of the universal entities within the vector that holds all the entities (owned, shared, and universal)
  _topologyToEntityContainerMap.find( p_topology )->second->SetEntityOwnershipIndexes( NSMesh::Ownership::o_UNIVERSAL, total_num_universal_entities, total_num_owned_entities + total_num_shared_entities );

  // Owned entities
  _populateEntityContainer( p_rank, _metaData->locally_owned_part(), p_topology, 0 );

  // Shared entities
  _populateEntityContainer( p_rank, _metaData->globally_shared_part() & !_metaData->locally_owned_part(), p_topology, total_num_owned_entities );

  // Universal entities
  _populateEntityContainer( p_rank, _metaData->universal_part() & (!_metaData->locally_owned_part()) & (!_metaData->globally_shared_part()), p_topology, total_num_owned_entities + total_num_shared_entities );
}


void
STKMesh::_populateEntityContainer( const stk::mesh::EntityRank& p_rank, const stk::mesh::Selector& p_selector, const NSMesh::Topology& p_topology, unsigned p_offset ) {
  // Get the entity container to be filled
  const auto& entity_container = _topologyToEntityContainerMap.find( p_topology )->second;

  // Get the container of the owned vertices only
  std::vector<stk::mesh::Bucket*> buckets = _bulkData->get_buckets( p_rank, p_selector );

  // Loop over all the entities
  size_t cnt = 0;
  for( size_t b=0; b<buckets.size(); ++b ) {
    stk::mesh::Bucket &bucket = *buckets[b];

    // For each entity
    for( size_t n=0; n<bucket.size(); ++n ) {
      stk::mesh::Entity stk_entity = bucket[n];

      // Create the NSMesh::NSEntity::Entity
      auto entity = std::make_unique<NSMesh::NSEntity::Entity>( p_topology, _bulkData->identifier( stk_entity ), p_offset + cnt );
      
      // Add the entity to the entity container
      entity_container->AddEntityAtIndex( entity, _bulkData->identifier( stk_entity ), p_offset + cnt );

      cnt++;
    }
  }
}


size_t
STKMesh::_getTotalNumEntities( const stk::mesh::EntityRank& p_rank, const stk::mesh::Selector& p_selector ) const {
  // Get the buckets containing the entities of a certain topology (i.e. rank) using a specified selector
  std::vector< stk::mesh::Bucket*> buckets = _bulkData->get_buckets( p_rank, p_selector );

  // Loop over all the buckets and get the total number of entities of a specified rank
  size_t result = 0;
  for( size_t b=0; b<buckets.size(); ++b ) {
    stk::mesh::Bucket &bucket = *buckets[b];
    result += bucket.size();
  }

  return result;
}


void 
STKMesh::_populateEntityConnectivityInformation() {
  // Ask the mesh query about the corresponding connectivity info to be stored
  const auto& connectivity_info_map = _meshQuery->GetConnectivityInfoProtocol();
  
  // Loop over all entity topologies
  for( const auto& it : connectivity_info_map ) {
    // Get the entity topology for which we are building the connectivity
    const NSMesh::Topology& from_entity = it.first;

    // For each entity topology, get all the entity topologies with which the connectivity relationship is explicitly stored
    const std::vector<NSMesh::Topology>& to_entity_vec = it.second;

    for( size_t i=0; i<to_entity_vec.size(); ++i ) {
      // For each entity with which an explicit connectivity information will be stored call the function that will store this connectivity relationship
      // For owned entities
      _populateEntityToEntityConnectivity( from_entity, to_entity_vec[i], _metaData->locally_owned_part() &(!_metaData->globally_shared_part()) );

      // For shared entities
      _populateEntityToEntityConnectivity( from_entity, to_entity_vec[i], _metaData->globally_shared_part() );

      // For universal entities
      _populateEntityToEntityConnectivity( from_entity, to_entity_vec[i], _metaData->universal_part() & (!_metaData->locally_owned_part()) & (!_metaData->globally_shared_part()) );
    }
  }
}


void
STKMesh::_populateEntityToEntityConnectivity( const NSMesh::Topology& p_fromEntityTopology, const NSMesh::Topology& p_toEntityTopology, const stk::mesh::Selector& p_selector ) {
  // Get the equivalent stk topologies of the entities we are building the connectivity for
  const stk::mesh::EntityRank& from_entity_topology_stk = _topologyToSTKMeshTopology.find( p_fromEntityTopology )->second;
  const stk::mesh::EntityRank& to_entity_topology_stk = _topologyToSTKMeshTopology.find( p_toEntityTopology )->second;

  // Get the container of the owned entities only in the stk data structure
  std::vector< stk::mesh::Bucket*> stk_buckets = _bulkData->get_buckets( from_entity_topology_stk, p_selector );

  // Loop over all the entities
  size_t cnt = 0;
  for( size_t b=0; b<stk_buckets.size(); ++b ) {
    stk::mesh::Bucket &stk_bucket = *stk_buckets[b];

    // For each stk entity
    for( size_t n=0; n<stk_bucket.size(); ++n ) {
      // Get the stk entity
      const stk::mesh::Entity& stk_entity = stk_bucket[n];

      // Get the entity
      NSMesh::NSEntity::Entity& from_entity = _getEntity( p_fromEntityTopology, _bulkData->identifier( stk_entity ) );

      // Get the stk connectivity
      const stk::mesh::Entity* entity_rel = _bulkData->begin( stk_entity, to_entity_topology_stk );
      
      for( size_t i=0; i<_bulkData->num_connectivity( stk_entity, to_entity_topology_stk ); ++i ) {
	// Get the stk entity with which to establish the connectivity relation
	const stk::mesh::Entity& to_entity_stk = entity_rel[i];

	// Get the entity with which to establish the connectivity relation
	const NSMesh::NSEntity::Entity& to_entity = GetEntity( p_toEntityTopology, _bulkData->identifier( to_entity_stk ) );

	// Add the connectivity relationship
	from_entity.AddConnectivity( p_toEntityTopology, &to_entity );
      }

      cnt++;
    }
  }
}


template<typename T>
void
STKMesh::_createSTKMeshParts( const T& p_nameToMeshPartMap ) {
  // Get the stk parts vector
  const stk::mesh::PartVector& stk_mesh_part_vec = _metaData->get_parts();
  
  // Loop over all the mesh parts
  for( const auto& it : p_nameToMeshPartMap ) {
    // For each mesh part
    const auto& mesh_part = it.second;

    // Get the mesh part topology
    const NSMesh::Topology& topology = mesh_part->GetPrimaryTopology();

    // Get the mesh part name
    const std::string& mesh_part_name = mesh_part->GetName();

    // Get the mesh part flag
    const int& id = mesh_part->GetId();

    // Check if this part has been defined in the mesh file and read by the stk mesh
    bool found = false;
    for( size_t p=0; p<stk_mesh_part_vec.size(); ++p ) {
      // Get the stk topology
      stk::mesh::EntityRank stk_topology = _topologyToSTKMeshTopology.find( topology )->second;
      
      if( ( id == stk_mesh_part_vec[p]->id() ) && ( stk_topology == stk_mesh_part_vec[p]->primary_entity_rank() ) ) {
	// Add the stk mesh part to the topology to stk part map
	_partNameToStkMeshPart.insert( std::make_pair( mesh_part_name, stk_mesh_part_vec[p] ) );
	found = true;

	// Get the topology of the cells that make up the part
	// Assumption is that the part has a uniform element topology
	stk::mesh::Selector selector = _metaData->locally_owned_part() & *stk_mesh_part_vec[ p ];
	std::vector<stk::mesh::Bucket*> buckets = _bulkData->get_buckets( stk_topology, selector );
	if( buckets.size() > 0 ) {
	  stk::mesh::Bucket& bucket = *buckets[0];
	  shards::CellTopology ct = stk::mesh::get_cell_topology( _bulkData->bucket( bucket[0] ).topology() );	  
	  _partNameToMeshCellTopologyInfo.insert( std::make_pair( it.first, std::make_pair( _getEltypeFromCellTopology( ct ), ct.getNodeCount() ) ) );
	}
	else {
	  // This means that the processor doesn't own any element of this part
	  // These default values used below will be overwritten with information received from other processors that own elements of the part
	  _partNameToMeshCellTopologyInfo.insert( std::make_pair( it.first, std::make_pair( NSMesh::ElementTopology::et_UNDEFINED, 0 ) ) );
	}
	break;
      }
    }
    
    // If the mesh part has not been defined in the mesh file through some flag for example throw an error
    NSExceptions::InvalidArgument( !found, "The mesh part " + mesh_part_name + " has not been defined in the mesh file " );
  }
}
    

void
STKMesh::_createSTKMeshParts( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  // Create the stk interior mesh parts
  _createSTKMeshParts( _nameToInteriorMeshPartMap );

  // Create the stk boundary mesh parts
  _createSTKMeshParts( _nameToBoundaryMeshPartMap );  

  // Parallel communicate the information about the element topology of each part to all the other processors
  _communicateMeshCellTopologyInfo( p_com );
}


void
STKMesh::_communicateMeshCellTopologyInfo( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  // Loop over the part name to cell info topology map
  for( auto& elem : _partNameToMeshCellTopologyInfo ) {
    // Get the part name
    std::string part_name = elem.first;

    // Get the element topology type
    NSMesh::ElementTopology element_topology = std::get<0>( elem.second );

    // Get the number of mesh vertices per element
    int num_mesh_vertices = std::get<1>( elem.second );
    
    // Communicate the element topology type a first time
    NSMesh::ElementTopology received_element_topology;
    p_com->AllReduce( &element_topology, &received_element_topology, 1, MPI_INT, MPI_MAX );
    
    // Send the number of vertices per mesh element type
    int received_num_mesh_vertices;
    p_com->AllReduce( &num_mesh_vertices, &received_num_mesh_vertices, 1, MPI_INT, MPI_MAX );

    // Set the element topology type and number of mesh vertices per element based on the received data
    NSMesh::ElementTopology& elem_topo = std::get<0>( elem.second );
    elem_topo = received_element_topology;
    
    int& num_vertices_per_mesh_element = std::get<1>( elem.second );
    num_vertices_per_mesh_element = received_num_mesh_vertices;
  }  
}


void
STKMesh::_populateMeshPartsEntities() {
  // Populate the interior mesh parts
  _populateInteriorMeshPartsEntities();

  // Populate the boundary mesh parts
  _populateBoundaryMeshPartsEntities();
}

    
void
STKMesh::_populateInteriorMeshPartsEntities() {
  // Loop over each interior mesh part
  for( const auto& it_mesh_part : _nameToInteriorMeshPartMap ) {
    // For each mesh part
    const std::unique_ptr<NSMesh::InteriorMeshPart>& mesh_part = it_mesh_part.second;

    // Get the mesh part name
    const std::string& mesh_part_name = it_mesh_part.first;

    // Get the stk part
    NSExceptions::LogicError( _partNameToStkMeshPart.find( mesh_part_name ) == _partNameToStkMeshPart.end(), "Trying to access an stk part that has not been defined" );
    const stk::mesh::Part& stk_mesh_part = *(_partNameToStkMeshPart.find( mesh_part_name )->second);

    // Set the topology information of the elements that make up the part, i.e. whether it is a tri, quad, tet, hex and how many mesh vertices per element
    const auto& element_topology = std::get<0>( _partNameToMeshCellTopologyInfo.find( it_mesh_part.first )->second );
    const auto& num_vertices_per_element = std::get<1>( _partNameToMeshCellTopologyInfo.find( it_mesh_part.first )->second );
    it_mesh_part.second->SetElementTopologyInfo( num_vertices_per_element, element_topology );

    // Create the vector holding the topology of the entities whose information will be stored in the interior mesh part
    std::vector<NSMesh::Topology> all_topology_vector = {NSMesh::Topology::tp_ELEMENT, NSMesh::Topology::tp_VERTEX, NSMesh::Topology::tp_EDGE};
    if( _ndim == 3 ) {
      all_topology_vector.push_back( NSMesh::Topology::tp_FACE );
    }
    else {
      // Store the data about the number of entities for this mesh part
      mesh_part->SetEntityStatsPerTopology( NSMesh::Topology::tp_FACE, NSMesh::Ownership::o_OWNED, 0 );
      mesh_part->SetEntityStatsPerTopology( NSMesh::Topology::tp_FACE, NSMesh::Ownership::o_SHARED, 0 );
      mesh_part->SetEntityStatsPerTopology( NSMesh::Topology::tp_FACE, NSMesh::Ownership::o_UNIVERSAL, 0 );
      
      // Get and set the total number of this entity on the all the processors
      mesh_part->SetGlobalNumEntityPerTopology( NSMesh::Topology::tp_FACE, 0 );
    }

    // Loop over each topology and store the information about entities of that topology in the interior mesh part
    for( size_t t=0; t<all_topology_vector.size(); ++t ) {
      _populateInteriorMeshPartsEntitiesPerTopology( mesh_part.get(), stk_mesh_part, all_topology_vector[t] );
    }
  }
}

    
void
STKMesh::_populateInteriorMeshPartsEntitiesPerTopology( NSMesh::InteriorMeshPart* const p_meshPart, const stk::mesh::Part& p_stkMeshPart,
							const NSMesh::Topology& p_topology ) {
  // Get the mesh part name
  auto mesh_part_name = p_meshPart->GetName();
  
  // Get the equivalent stk topology
  stk::mesh::EntityRank stk_topology = _topologyToSTKMeshTopology.find( p_topology )->second;
  
  // Get the total number of entities of a specified topology defined on the mesh part
  size_t total_num_entities_part = _getTotalNumEntities( stk_topology, p_stkMeshPart );
  
  // Allocate memory for the part
  // This call below returns all the entities of a specified topology (owned, shared, universal)
  p_meshPart->AllocateMemoryPerTopology( p_topology, total_num_entities_part );
  
  // Get the total number of owned, shared, and universal entities per mesh part
  stk::mesh::Selector owned_selector = _metaData->locally_owned_part() & p_stkMeshPart;
  size_t total_num_owned_entities_part = _getTotalNumEntities( stk_topology, owned_selector );
  stk::mesh::Selector shared_selector = (_metaData->globally_shared_part() & !_metaData->locally_owned_part()) & p_stkMeshPart;
  size_t total_num_shared_entities_part = _getTotalNumEntities( stk_topology, shared_selector );
  stk::mesh::Selector universal_selector = (_metaData->universal_part() & (!_metaData->locally_owned_part()) & (!_metaData->globally_shared_part())) & p_stkMeshPart;
  size_t total_num_universal_entities_part = _getTotalNumEntities( stk_topology, universal_selector );
  
  // Store the data about the number of entities for this mesh part
  p_meshPart->SetEntityStatsPerTopology( p_topology, NSMesh::Ownership::o_OWNED, total_num_owned_entities_part );
  p_meshPart->SetEntityStatsPerTopology( p_topology, NSMesh::Ownership::o_SHARED, total_num_shared_entities_part );
  p_meshPart->SetEntityStatsPerTopology( p_topology, NSMesh::Ownership::o_UNIVERSAL, total_num_universal_entities_part );
  
  // Get and set the total number of this entity, defined on the mesh part, on the all the processors
  size_t global_num_entities_part = 0;
  _com->AllReduce( &total_num_owned_entities_part, &global_num_entities_part, 1, my_MPI_SIZE_T, MPI_SUM );
  p_meshPart->SetGlobalNumEntityPerTopology( p_topology, global_num_entities_part );
  
  // Set the starting and ending index of the owned entities within the vector that holds all the entities (owned, shared, and universal)
  p_meshPart->SetEntityOwnershipIndexesPerTopology( p_topology, NSMesh::Ownership::o_OWNED, total_num_owned_entities_part, 0 );
  
  // Set the starting and ending index of the shared entities within the vector that holds all the entities (owned, shared, and universal)
  p_meshPart->SetEntityOwnershipIndexesPerTopology( p_topology, NSMesh::Ownership::o_SHARED, total_num_shared_entities_part, total_num_owned_entities_part );
  
  // Set the starting and ending index of the universal entities within the vector that holds all the entities (owned, shared, and universal)
  p_meshPart->SetEntityOwnershipIndexesPerTopology( p_topology, NSMesh::Ownership::o_UNIVERSAL, total_num_universal_entities_part, total_num_owned_entities_part + total_num_shared_entities_part );
  
  // Populate the owned entities gids vector of the mesh part
  _populateInteriorMeshPartEntityVector( mesh_part_name, stk_topology, owned_selector, p_topology, 0 );
  
  // Populate the shared entities gids vector of the mesh part
  _populateInteriorMeshPartEntityVector( mesh_part_name, stk_topology, shared_selector, p_topology, total_num_owned_entities_part );
  
  // Populate the univeral entities gids vector of the mesh part
  _populateInteriorMeshPartEntityVector( mesh_part_name, stk_topology, universal_selector, p_topology, total_num_owned_entities_part + total_num_shared_entities_part );
}


void
STKMesh::_populateInteriorMeshPartEntityVector( const std::string& p_name, const stk::mesh::EntityRank& p_stkTopology, const stk::mesh::Selector& p_selector, const NSMesh::Topology& p_topology, unsigned p_offset ) {
  // Get the mesh part whose entities vector will be filled
  auto& meshPart = *_nameToInteriorMeshPartMap.find( p_name )->second;

  // Get the buckets containing the entity of a particular topology and based on the selector
  std::vector<stk::mesh::Bucket*> buckets = _bulkData->get_buckets( p_stkTopology, p_selector );
  
  // Loop over all the entities
  size_t cnt = 0;
  for( size_t b=0; b<buckets.size(); ++b ) {
    stk::mesh::Bucket &bucket = *buckets[b];

    // For each entity
    for( size_t n=0; n<bucket.size(); ++n ) {
      // Get the stk entity
      const stk::mesh::Entity& stk_entity = bucket[n];

      // Get the pointer to the NSMesh::NSEntity::Entity
      // This object has been defined before the call to this function
      auto& entity = _getEntity( p_topology, _bulkData->identifier( stk_entity ) );
      
      // Add the entity lid to the mesh part lids vector
      meshPart.AddEntityAtIndex( &entity, p_topology, _bulkData->identifier( stk_entity ), p_offset + cnt );

      cnt++;
    }
  }
}


void
STKMesh::_populateBoundaryMeshPartsEntities() {
  // Loop over each boundary mesh part
  for( const auto& it_mesh_part : _nameToBoundaryMeshPartMap ) {
    // For each mesh part
    const std::unique_ptr<NSMesh::BoundaryMeshPart>& mesh_part = it_mesh_part.second;

    // Get the mesh part name
    const std::string& mesh_part_name = it_mesh_part.first;

    // Get the stk part
    NSExceptions::LogicError( _partNameToStkMeshPart.find( mesh_part_name ) == _partNameToStkMeshPart.end(), "Trying to access an stk part that has not been defined" );
    const stk::mesh::Part& stk_mesh_part = *(_partNameToStkMeshPart.find( mesh_part_name )->second);

    // Set the topology information of the elements that make up the part, i.e. whether it is a tri, quad, tet, hex and how many mesh vertices per element
    const auto& element_topology = std::get<0>( _partNameToMeshCellTopologyInfo.find( it_mesh_part.first )->second );
    const auto& num_vertices_per_element = std::get<1>( _partNameToMeshCellTopologyInfo.find( it_mesh_part.first )->second );
    it_mesh_part.second->SetElementTopologyInfo( num_vertices_per_element, element_topology );
    
    // Populate the boundary mesh part
    _populateBoundaryMeshPartsEntitiesForAllTopologies( mesh_part.get(), stk_mesh_part );
  }
}


void
STKMesh::_populateBoundaryMeshPartsEntitiesForAllTopologies( NSMesh::BoundaryMeshPart* const p_meshPart, const stk::mesh::Part& p_stkMeshPart ) {
  // Get the mesh part name
  auto mesh_part_name = p_meshPart->GetName();
  
  // Get the equivalent stk topology
  stk::mesh::EntityRank stk_topology = _topologyToSTKMeshTopology.find( GetClosureTopology() )->second;
  
  // Get the selectors corresponding to the different ownership types
  stk::mesh::Selector owned_selector = _metaData->locally_owned_part() & p_stkMeshPart;
  stk::mesh::Selector shared_selector = (_metaData->globally_shared_part() & !_metaData->locally_owned_part()) & p_stkMeshPart;
  stk::mesh::Selector universal_selector = (_metaData->universal_part() & (!_metaData->locally_owned_part()) & (!_metaData->globally_shared_part())) & p_stkMeshPart;
  
  // Populate the owned entities gids vector of the mesh part
  std::map<std::string, std::vector<unsigned> > owned_interior_part_name_to_closure_entity_gid;
  std::map<std::string, std::vector<unsigned> > owned_interior_part_name_to_element_gid;
  _populateBoundaryMeshPartEntityVector( stk_topology, owned_selector, owned_interior_part_name_to_closure_entity_gid, owned_interior_part_name_to_element_gid );
  
  // Populate the shared entities gids vector of the mesh part
  std::map<std::string, std::vector<unsigned> > shared_interior_part_name_to_closure_entity_gid;
  std::map<std::string, std::vector<unsigned> > shared_interior_part_name_to_element_gid;
  _populateBoundaryMeshPartEntityVector( stk_topology, shared_selector, shared_interior_part_name_to_closure_entity_gid, shared_interior_part_name_to_element_gid );
  
  // Populate the univeral entities gids vector of the mesh part
  std::map<std::string, std::vector<unsigned> > universal_interior_part_name_to_closure_entity_gid;
  std::map<std::string, std::vector<unsigned> > universal_interior_part_name_to_element_gid;
  _populateBoundaryMeshPartEntityVector( stk_topology, universal_selector, universal_interior_part_name_to_closure_entity_gid, universal_interior_part_name_to_element_gid );

  // Make sure that all three maps above have the same size on each processor
  _synchronizeBoundaryMeshMapsPerTopologyPerProcessor( owned_interior_part_name_to_closure_entity_gid,
						       shared_interior_part_name_to_closure_entity_gid,
						       universal_interior_part_name_to_closure_entity_gid );
  _synchronizeBoundaryMeshMapsPerTopologyPerProcessor( owned_interior_part_name_to_element_gid,
						       shared_interior_part_name_to_element_gid,
						       universal_interior_part_name_to_element_gid );

  // Make sure that all three maps above have the same size across processors
  _synchronizeBoundaryMeshMapsPerTopologyPerProcessor( owned_interior_part_name_to_closure_entity_gid,
						       shared_interior_part_name_to_closure_entity_gid,
						       universal_interior_part_name_to_closure_entity_gid,
						       owned_interior_part_name_to_element_gid,
						       shared_interior_part_name_to_element_gid,
						       universal_interior_part_name_to_element_gid );
  
  // Get the boundary mesh part whose entities vector will be filled
  auto& mesh_part = *_nameToBoundaryMeshPartMap.find( mesh_part_name )->second;
  
  // Loop over each interior part attached to the boundary part
  std::vector<std::string> attached_interior_mesh_parts_name;
  for( const auto& interior_part : owned_interior_part_name_to_closure_entity_gid ) {
    // Get the interior part name
    const auto& interior_part_name = interior_part.first;

    // Get the vectors holding the owned, shared, and universal closure entity gids
    const auto& owned_closure_entity_gid = owned_interior_part_name_to_closure_entity_gid.find( interior_part_name )->second;
    const auto& shared_closure_entity_gid = shared_interior_part_name_to_closure_entity_gid.find( interior_part_name )->second;
    const auto& universal_closure_entity_gid = universal_interior_part_name_to_closure_entity_gid.find( interior_part_name )->second;

    // Complete the storing of the data about the closure entities for the current interior mesh part
    _completeBoundaryMeshPartPerTopology( interior_part_name, mesh_part, owned_closure_entity_gid,
					  shared_closure_entity_gid, universal_closure_entity_gid,
					  GetClosureTopology() );

    // Get the vectors holding the owned, shared, and universal element gids
    const auto& owned_element_gid = owned_interior_part_name_to_element_gid.find( interior_part_name )->second;
    const auto& shared_element_gid = shared_interior_part_name_to_element_gid.find( interior_part_name )->second;
    const auto& universal_element_gid = universal_interior_part_name_to_element_gid.find( interior_part_name )->second;

    // Complete the storing of the data about the elemnts for the current interior mesh part
    _completeBoundaryMeshPartPerTopology( interior_part_name, mesh_part, owned_element_gid,
					  shared_element_gid, universal_element_gid,
					  NSMesh::Topology::tp_ELEMENT );

    // Check if the interior part is attached to the boundary part
    if( owned_closure_entity_gid.size() > 0 || shared_closure_entity_gid.size() > 0 || universal_closure_entity_gid.size() > 0 ||
	owned_element_gid.size() > 0 || shared_element_gid.size() > 0 || universal_element_gid.size() > 0 ) {
      attached_interior_mesh_parts_name.push_back( interior_part_name );
    }
  }

  // Set the list of interior part names attached to this boundary mesh part
  p_meshPart->SetAttachedInteriorMeshPartsName( attached_interior_mesh_parts_name );
}


void
STKMesh::_synchronizeBoundaryMeshMapsPerTopologyPerProcessor( std::map<std::string, std::vector<unsigned> >& p_ownedData,
							      std::map<std::string, std::vector<unsigned> >& p_sharedData,
							      std::map<std::string, std::vector<unsigned> >& p_universalData ) const {
  // Synchrnize the maps locally
  for( const auto& gids_data : p_ownedData ) {
    if( p_sharedData.find( gids_data.first ) == p_sharedData.end() ) {
      p_sharedData.insert( std::make_pair( gids_data.first, std::vector<unsigned>() ) );
    }
    if( p_universalData.find( gids_data.first ) == p_universalData.end() ) {
      p_universalData.insert( std::make_pair( gids_data.first, std::vector<unsigned>() ) );
    }
  }

  for( const auto& gids_data : p_sharedData ) {
    if( p_ownedData.find( gids_data.first ) == p_ownedData.end() ) {
      p_ownedData.insert( std::make_pair( gids_data.first, std::vector<unsigned>() ) );
    }
    if( p_universalData.find( gids_data.first ) == p_universalData.end() ) {
      p_universalData.insert( std::make_pair( gids_data.first, std::vector<unsigned>() ) );
    }
  }

  for( const auto& gids_data : p_universalData ) {
    if( p_ownedData.find( gids_data.first ) == p_ownedData.end() ) {
      p_ownedData.insert( std::make_pair( gids_data.first, std::vector<unsigned>() ) );
    }
    if( p_sharedData.find( gids_data.first ) == p_sharedData.end() ) {
      p_sharedData.insert( std::make_pair( gids_data.first, std::vector<unsigned>() ) );
    }
  }
}


void
STKMesh::_synchronizeBoundaryMeshMapsPerTopologyPerProcessor( std::map<std::string, std::vector<unsigned> >& p_closureEntityOwnedData,
							      std::map<std::string, std::vector<unsigned> >& p_closureEntitySharedData,
							      std::map<std::string, std::vector<unsigned> >& p_closureEntityUniversalData,
							      std::map<std::string, std::vector<unsigned> >& p_elementOwnedData,
							      std::map<std::string, std::vector<unsigned> >& p_elementSharedData,
							      std::map<std::string, std::vector<unsigned> >& p_elementUniversalData ) const {
  // If we have multiple processes
  if( _com->GetSize() > 1 ) {
    auto com_size = _com->GetSize();
    std::vector<int> num_to_receive( com_size, 0 );

    // Determine how many interior part names each processor will send to the other processes
    num_to_receive[ _com->GetRank() ] = static_cast<int>( p_closureEntityOwnedData.size() );
    for( int i=0; i<com_size; ++i ) {
      _com->Broadcast( &num_to_receive[ i ], 1, MPI_INT, i );
    }

    // Determine the sizes of the strings to be communicated by each processor
    std::vector<std::vector<int> > to_receive_string_pos( com_size );
    for( int i=0; i<static_cast<int>( num_to_receive.size() ); ++i ) {
      if( i == _com->GetRank() ) {
	for( const auto& string_data : p_closureEntityOwnedData ) {
	  int position = std::distance( _nameToInteriorMeshPartMap.begin(), _nameToInteriorMeshPartMap.find( string_data.first ) );
	  to_receive_string_pos[ i ].push_back( position );
	}
      }
      else {
	if( num_to_receive[ i ] > 0 ) {
	  to_receive_string_pos[ i ].resize( num_to_receive[ i ] );
	}
      }
    }

    for( size_t i=0; i<to_receive_string_pos.size(); ++i ) {
      _com->Broadcast( &to_receive_string_pos[i][0], to_receive_string_pos[i].size(), MPI_INT, i );
    }
    
    // Update the owned, shared and universal data maps with the communicated data if necessary
    for( size_t i=0; i<to_receive_string_pos.size(); ++i ) {
      for( size_t j=0; j<to_receive_string_pos[i].size(); ++j ) {
	// Get the interior part name
	auto it = _nameToInteriorMeshPartMap.begin();
	std::advance( it, to_receive_string_pos[i][j] );
	auto interior_part_name = it->first;

	// Synchronize the closure entity data maps
	if( p_closureEntityOwnedData.find( interior_part_name ) == p_closureEntityOwnedData.end() ) {
	  p_closureEntityOwnedData.insert( std::make_pair( interior_part_name, std::vector<unsigned>() ) );
	}	
	if( p_closureEntitySharedData.find( interior_part_name ) == p_closureEntitySharedData.end() ) {
	  p_closureEntitySharedData.insert( std::make_pair( interior_part_name, std::vector<unsigned>() ) );
	}	
	if( p_closureEntityUniversalData.find( interior_part_name ) == p_closureEntityUniversalData.end() ) {
	  p_closureEntityUniversalData.insert( std::make_pair( interior_part_name, std::vector<unsigned>() ) );
	}

	// Synchronize the element data maps
	if( p_elementOwnedData.find( interior_part_name ) == p_elementOwnedData.end() ) {
	  p_elementOwnedData.insert( std::make_pair( interior_part_name, std::vector<unsigned>() ) );
	}	
	if( p_elementSharedData.find( interior_part_name ) == p_elementSharedData.end() ) {
	  p_elementSharedData.insert( std::make_pair( interior_part_name, std::vector<unsigned>() ) );
	}	
	if( p_elementUniversalData.find( interior_part_name ) == p_elementUniversalData.end() ) {
	  p_elementUniversalData.insert( std::make_pair( interior_part_name, std::vector<unsigned>() ) );
	}
      }
    }
  }
}
    
    
void
STKMesh::_populateBoundaryMeshPartEntityVector( const stk::mesh::EntityRank& p_stkTopology, const stk::mesh::Selector& p_selector,
						std::map<std::string, std::vector<unsigned> >& p_interiorPartNameToClosureEntityGid,
						std::map<std::string, std::vector<unsigned> >& p_interiorPartNameToElementGid ) {
  // Get the buckets containing the closure entities and based on the selector
  std::vector<stk::mesh::Bucket*> buckets = _bulkData->get_buckets( p_stkTopology, p_selector );
  
  // Loop over all the closure entities
  for( size_t b=0; b<buckets.size(); ++b ) {
    stk::mesh::Bucket &bucket = *buckets[b];

    // For each closure entity
    for( size_t n=0; n<bucket.size(); ++n ) {
      // Get the stk closure entity
      const stk::mesh::Entity& stk_closure_entity = bucket[n];

      // Get the gid of the closure entity
      unsigned closure_entity_gid = _bulkData->identifier( stk_closure_entity );

      // Get the connected element to this closure entity
      const stk::mesh::Entity* stk_element = _bulkData->begin( stk_closure_entity, stk::topology::ELEMENT_RANK );

      // Check that only one element is connected to the closure entity
      NSExceptions::LogicError( _bulkData->num_connectivity( stk_closure_entity, stk::topology::ELEMENT_RANK ) != 1, "Only one element can be connected to a closure entity" );
      
      // Get the gid of the connected element
      unsigned stk_element_gid = _bulkData->identifier( *stk_element );

      // Get the id of the interior mesh part to which the element belongs
      int interior_mesh_part_id = -1;
      for( const auto& interior_mesh_part : _nameToInteriorMeshPartMap ) {
	if( interior_mesh_part.second->IsPresent( stk_element_gid, NSMesh::Topology::tp_ELEMENT ) ) {
	  // Set the interior mesh part id
	  interior_mesh_part_id = interior_mesh_part.second->GetId();

	  // Add the data to the corresponding maps
	  if( p_interiorPartNameToElementGid.find( interior_mesh_part.first ) != p_interiorPartNameToElementGid.end() ) {
	    p_interiorPartNameToElementGid.find( interior_mesh_part.first )->second.push_back( stk_element_gid );
	    p_interiorPartNameToClosureEntityGid.find( interior_mesh_part.first )->second.push_back( closure_entity_gid );
	  }
	  else {
	    p_interiorPartNameToElementGid.insert( std::make_pair( interior_mesh_part.first, std::vector<unsigned>( 1, stk_element_gid ) ) );
	    p_interiorPartNameToClosureEntityGid.insert( std::make_pair( interior_mesh_part.first, std::vector<unsigned>( 1, closure_entity_gid ) ) );;
	  }
	  break;
	}
      }

      // Check that the element belongs to an interior mesh part
      NSExceptions::LogicError( interior_mesh_part_id == -1, "The element attached to the boundary mesh part must belong to an interior mesh part" );
    }
  }
}


void
STKMesh::_completeBoundaryMeshPartPerTopology( const std::string& p_interiorPartName, NSMesh::BoundaryMeshPart& p_boundaryMeshPart,
					       const std::vector<unsigned>& p_ownedEntityGid, const std::vector<unsigned>& p_sharedEntityGid,
					       const std::vector<unsigned>& p_universalEntityGid, const NSMesh::Topology& p_topology ) {
  // Allocate memory for the entity attached to the interior mesh part
  unsigned total_number_entities = p_ownedEntityGid.size() + p_sharedEntityGid.size() + p_universalEntityGid.size();
  p_boundaryMeshPart.AllocateMemoryPerTopology( p_topology, p_interiorPartName, total_number_entities );

  // Store the data about the number of entities for the interior part attached to the boundary part
  p_boundaryMeshPart.SetEntityStatsPerTopology( p_topology, p_interiorPartName, NSMesh::Ownership::o_OWNED, p_ownedEntityGid.size() );
  p_boundaryMeshPart.SetEntityStatsPerTopology( p_topology, p_interiorPartName, NSMesh::Ownership::o_SHARED, p_sharedEntityGid.size() );
  p_boundaryMeshPart.SetEntityStatsPerTopology( p_topology, p_interiorPartName, NSMesh::Ownership::o_UNIVERSAL, p_universalEntityGid.size() );

  // Get and set the total number of this entity, defined on the mesh part, on the all the processors
  size_t global_num_entities_part = 0;
  _com->AllReduce( &total_number_entities, &global_num_entities_part, 1, my_MPI_SIZE_T, MPI_SUM );
  p_boundaryMeshPart.SetGlobalNumEntityPerTopology( p_topology, p_interiorPartName, global_num_entities_part );
  
  // Set the starting and ending index of the owned entities within the vector that holds all the entities (owned, shared, and universal)
  p_boundaryMeshPart.SetEntityOwnershipIndexesPerTopology( p_topology, p_interiorPartName, NSMesh::Ownership::o_OWNED,
							   p_ownedEntityGid.size(), 0 );
  _addEntityGidToBoundaryMeshPart( p_ownedEntityGid, p_boundaryMeshPart, p_interiorPartName, 0, p_topology );

  // Set the starting and ending index of the shared entities within the vector that holds all the entities (owned, shared, and universal)
  p_boundaryMeshPart.SetEntityOwnershipIndexesPerTopology( p_topology, p_interiorPartName, NSMesh::Ownership::o_SHARED, p_sharedEntityGid.size(),
							   p_ownedEntityGid.size() );
  _addEntityGidToBoundaryMeshPart( p_sharedEntityGid, p_boundaryMeshPart, p_interiorPartName, p_ownedEntityGid.size(), p_topology );

  // Set the starting and ending index of the universal entities within the vector that holds all the entities (owned, shared, and universal)
  p_boundaryMeshPart.SetEntityOwnershipIndexesPerTopology( p_topology, p_interiorPartName, NSMesh::Ownership::o_UNIVERSAL, p_universalEntityGid.size(),
							   p_ownedEntityGid.size() + p_sharedEntityGid.size() );
  _addEntityGidToBoundaryMeshPart( p_universalEntityGid, p_boundaryMeshPart, p_interiorPartName, p_ownedEntityGid.size() + p_sharedEntityGid.size(), p_topology );
}


void
STKMesh::_addEntityGidToBoundaryMeshPart( const std::vector<unsigned>& p_entityGids, NSMesh::BoundaryMeshPart& p_meshPart,
					  const std::string& p_interiorPartName, unsigned p_offset, const NSMesh::Topology& p_topology ) {
  // Loop over the entity gids vector
  size_t cnt = 0;
  for( const auto& entity_gid : p_entityGids ) {
    // Get the entity
    auto& entity = _getEntity( p_topology, entity_gid );
    
    // Add the entity lid to the mesh part lids vector
    p_meshPart.AddEntityAtIndex( &entity, p_interiorPartName, p_topology, entity_gid, p_offset + cnt );
    cnt++;
  }
}

    
void
STKMesh::PopulateCoordinatesField( const NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*>& p_bucket, const Kokkos::View<double***>& p_coordinates ) const {
  // Get the stk coordinates vector
  stk::mesh::FieldBase const* stkCoordinatesField = _metaData->coordinate_field();
  
  // Loop over the elements
  for( size_t e=0; e<p_bucket.GetSize(); ++e ) {
    // Get the element
    const auto& element = *p_bucket[e];
    
    // Get the vertices attached to the element
    std::vector<const NSMesh::NSEntity::Entity*> attached_vertices;
    Begin( NSMesh::Topology::tp_VERTEX, element, attached_vertices );
    
    // Loop over the vertices
    for( size_t n=0; n<attached_vertices.size(); ++n ) {
      // Get the gid of the vertex before the renumbering
      const auto& old_vertex_gid = _newToOldVerticesGIDMap.find( attached_vertices[n]->GetGid() )->second;
      
      // Get the stk vertex
      const stk::mesh::Entity& vertex = _bulkData->get_entity( stk::topology::NODE_RANK, old_vertex_gid );
      
      // Get the vertex coordinates from the stk field
      double* stkValues = (double*)stk::mesh::field_data( *stkCoordinatesField, vertex );

      // Set the values of the mesh coordinates field from the stk coordinates field
      for( int dof=0; dof<_ndim; ++dof ) {
	p_coordinates( e, n, dof ) = stkValues[dof];
      }
    }
  }
}


void
STKMesh::_populateEntityComInfo() {
  // Get entities living at the interface between the different processors
  std::vector<stk::mesh::Bucket*> buckets = _bulkData->get_buckets( stk::topology::NODE_RANK, _metaData->globally_shared_part() );
  
  // Loop over the buckets
  std::vector<stk::mesh::Entity> allVertices;
  std::vector<stk::mesh::Entity> allElements;
  std::vector<stk::mesh::Entity> allEdges;
  std::vector<stk::mesh::Entity> allFaces;
  for( size_t b=0; b<buckets.size(); ++b ) {
    stk::mesh::Bucket& bucket = *buckets[b];
    
    // Loop over the vertices
    for( size_t n=0; n<bucket.size(); ++n ) {
      // Get the stk vertex
      stk::mesh::Entity vertex = bucket[n];

      // Get the elements attached
      const stk::mesh::Entity* elem_rels = _bulkData->begin_elements( vertex );

      // Loop over the elements
      for( size_t e=0; e<_bulkData->num_elements( vertex ); ++e ) {
	const stk::mesh::Entity& element = elem_rels[e];

	// Push back the elements
	allElements.push_back( element );

	// Get the vertices attached to the element
	const stk::mesh::Entity* node_rels = _bulkData->begin_nodes( element );

	// Push back the vertices to the all vertices vector
	for( size_t no=0; no<_bulkData->num_nodes( element ); ++no ) {
	  const stk::mesh::Entity& subnode = node_rels[no];
	  allVertices.push_back( subnode );
	}

	// Get the edges attached to the element
	const stk::mesh::Entity* edge_rels = _bulkData->begin_edges( element );
	
	// Push back the edges to the all edges vector
	for( size_t ed=0; ed<_bulkData->num_edges( element ); ++ed ) {
	  const stk::mesh::Entity& subedge = edge_rels[ed];
	  allEdges.push_back( subedge );
	}

	if( _ndim == 3 ) {
	  // Get the faces attached to the element
	  const stk::mesh::Entity* face_rels = _bulkData->begin_faces( element );
	  
	  // Push back the faces to the all faces vector
	  for( size_t fa=0; fa<_bulkData->num_faces( element ); ++fa ) {
	    const stk::mesh::Entity& subface = face_rels[fa];
	    allFaces.push_back( subface );
	  }
	}
      }
    }
  }

  // Remove the duplicated entities
  std::sort( allVertices.begin(), allVertices.end() );
  allVertices.erase( std::unique( allVertices.begin(), allVertices.end() ), allVertices.end() );

  // Populate vertex com info
  _populateEntityComInfoForRank( NSMesh::Topology::tp_VERTEX, allVertices );

  std::sort( allElements.begin(), allElements.end() );
  allElements.erase( std::unique( allElements.begin(), allElements.end() ), allElements.end() );

  // Populate element com info
  _populateEntityComInfoForRank( NSMesh::Topology::tp_ELEMENT, allElements );

  std::sort( allEdges.begin(), allEdges.end() );
  allEdges.erase( std::unique( allEdges.begin(), allEdges.end() ), allEdges.end() );

  // Populate edge com info
  _populateEntityComInfoForRank( NSMesh::Topology::tp_EDGE, allEdges );

  if( _ndim == 3 ) {
    std::sort( allFaces.begin(), allFaces.end() );
    allFaces.erase( std::unique( allFaces.begin(), allFaces.end() ), allFaces.end() );
    
    // Populate face com info
    _populateEntityComInfoForRank( NSMesh::Topology::tp_FACE, allFaces );
  }
}


void
STKMesh::_populateEntityComInfoForRank( const NSMesh::Topology& p_topology, const std::vector<stk::mesh::Entity>& p_entities ) {
  // Get the map to be filled
  _topoToEntityGidProcMap.insert( std::make_pair( p_topology, std::map<unsigned, NSMesh::NSFields::EntityComInfo>{} ) );
  auto& entityToProcsRankMap = _topoToEntityGidProcMap.find( p_topology )->second;

  // Get the rank of the processor
  int my_rank = _com->GetRank();
  
  // Loop over the entities
  for( size_t n=0; n<p_entities.size(); ++n ) {
    // Get the stk entity
    const stk::mesh::Entity& entity = p_entities[n];

    // Get the stk entity gid
    unsigned gid = _bulkData->identifier( entity );

    // Allocate memory for the entity com info object. This might not be needed eventually
    auto entity_com_info = NSMesh::NSFields::EntityComInfo();

    // Get the owning processor
    int owningProc = _bulkData->parallel_owner_rank( entity );
    entity_com_info.SetOwningProcID( owningProc );
    
    if( owningProc == my_rank ) {
      // Get the rank of the processors for which this entity is shared
      std::vector<int> sharingProcs;
      _bulkData->comm_shared_procs( _bulkData->entity_key( entity ), sharingProcs );
      
      // Get the rank of the processors for which this entity is universal
      std::vector<int> universalProcs;
      _bulkData->comm_procs( entity, universalProcs );
    
      // Remove the owning processor rank from the sharing procs vector
      std::vector<int>::iterator it_o = std::find( sharingProcs.begin(), sharingProcs.end(), owningProc );
      if( it_o != sharingProcs.end() ) {
	sharingProcs.erase( it_o );
      }
      
      // Remove the sharing processors from the vector of universal processors
      std::vector<int> indexToRemove;
      for( size_t p=0; p<universalProcs.size(); ++p ) {
	if( std::find( sharingProcs.begin(), sharingProcs.end(), universalProcs[p] ) != sharingProcs.end() ) {
	  indexToRemove.push_back( p );
	}
      }
      
      for( int p=indexToRemove.size() - 1; p>=0; --p ) {
	universalProcs.erase( universalProcs.begin() + indexToRemove[p] );
      }

      // Remove the owning processor rank from the universal procs vector
      it_o = std::find( universalProcs.begin(), universalProcs.end(), owningProc );
      if( it_o != universalProcs.end() ) {
	universalProcs.erase( it_o );
      }
    
      // Allocate memory for this entity in the corresponding map if needed
      if( sharingProcs.size() > 0 || universalProcs.size() > 0 ) {
	entity_com_info.SetSharingProcID( sharingProcs );
	entity_com_info.SetUniversalProcID( universalProcs );
	entityToProcsRankMap.insert( std::make_pair( gid, entity_com_info ) );
      }
    }
    else {
      // Get the rank of the processors for which this entity is shared
      std::vector<int> tmp_shared_procs_vec;
      _bulkData->comm_shared_procs( _bulkData->entity_key( entity ), tmp_shared_procs_vec );
      if( std::find( tmp_shared_procs_vec.begin(), tmp_shared_procs_vec.end(), owningProc ) != tmp_shared_procs_vec.end() ) {
	entity_com_info.AddSharingProcID( my_rank );
	entityToProcsRankMap.insert( std::make_pair( gid, entity_com_info ) );
      }
      else {
	std::vector<int> tmp_universal_procs_vec;
	_bulkData->comm_procs( entity, tmp_universal_procs_vec );
	if( std::find( tmp_universal_procs_vec.begin(), tmp_universal_procs_vec.end(), owningProc ) != tmp_universal_procs_vec.end() ) {
	  entity_com_info.AddUniversalProcID( my_rank );
	  entityToProcsRankMap.insert( std::make_pair( gid, entity_com_info ) );
	}
      }
    }
  }
}


void
STKMesh::ClearTPLMesh() {
  _meshReader.reset();

  _meshReader = nullptr;
  _bulkData = nullptr;
  _metaData = nullptr;
}


void
STKMesh::_renumberEntitiesImpl( const NSMesh::Topology& p_topology, const stk::mesh::EntityRank& p_rank, std::map<unsigned, unsigned>& p_oldToNewGIDMap, size_t p_offset ) {
  // Get the entity gid to processor map for the given topology
  // This data structure holds data about the entities defined at the interface between the processors
  // Those entities could be shared by two processors or owned by a processor and within the aura of another processor
  auto& gidToProcessorMap = _topoToEntityGidProcMap.find( p_topology )->second;

  // Get the number of owned entities on the processor
  size_t numOfOwnedEntities = _topologyToEntityContainerMap.find( p_topology )->second->GetNumberOfEntity( {NSMesh::Ownership::o_OWNED } );
  size_t startRange = 0;

  // Communicate the number of owned entities with the other processors
  // This step is important so that each processor will know what integer it must start from in generating the new numbers for the entities
  _com->ExclusiveScan( &numOfOwnedEntities, &startRange, 1, my_MPI_SIZE_T, MPI_SUM );

  // Offset the start range by the global number of entities of different topologies
  startRange += p_offset;
  
  // Get the bucket of owned entities
  std::vector<stk::mesh::Bucket*> buckets = _bulkData->get_buckets( p_rank, _metaData->locally_owned_part() );
  
  // Loop over the vector of buckets
  size_t cnt = startRange + 1;
  for( size_t b=0; b<buckets.size(); ++b ) {
    // For each bucket
    const stk::mesh::Bucket& bucket = *buckets[b];
    
    // Loop over the entities of each bucket
    for( size_t e=0; e<bucket.size(); ++e ) {
      // Get the stk entity
      const stk::mesh::Entity& entity = bucket[e];
      
      // Get a non-const version of the entity
      NSMesh::NSEntity::Entity& non_const_entity = _getEntity( p_topology, _bulkData->identifier( entity ) );

      // If the entity is defined at the interface between processors meaning that this entity is defined on multiple processors
      // Then store the old and new gids for this entity
      if( gidToProcessorMap.find( non_const_entity.GetGid() ) != gidToProcessorMap.end() ) {
	p_oldToNewGIDMap.insert( std::make_pair( non_const_entity.GetGid(), cnt ) );
      }

      // Fill the _newToOldVerticesGIDMap. This map is needed later to create a finite element field that holds the coordinates of the mesh vertices
      if( p_topology == NSMesh::Topology::tp_VERTEX ) {
	_newToOldVerticesGIDMap.insert( std::make_pair( cnt, non_const_entity.GetGid() ) );
      }
      
      // Set the new gid of the entity
      non_const_entity.SetGid( cnt );

      ++cnt;
    }
  }
}


void
STKMesh::_communicateNewGIDs( const NSMesh::Topology& p_topology, const std::map<unsigned, unsigned>& p_oldToNewGIDMap ) {
  // Check that this is a parallel simulation
  if( _com->IsParallel() ) {
    // Get the communicator size
    int com_size = _com->GetSize();
    
    // Get the processor rank
    int proc_rank = _com->GetRank();
    
    // Vector specifyng the number of data to be sent and received by each processor
    std::vector<unsigned> send_size( com_size, 0 );
    std::vector<unsigned> recv_size( com_size, 0 );
    std::vector<std::vector<unsigned> > to_send_gid( com_size, std::vector<unsigned>{} );
    
    // Get the map holding data about entities that are shared or universal for the different processors
    // At this stage, this map holds the old entity gid even on the processor on which this entity gid has been updated
    const auto& entity_com_info_map = _topoToEntityGidProcMap.find( p_topology )->second;
    
    // Loop over the entities
    for( const auto& it_e : entity_com_info_map ) {
      // Get the data structure holding the communication info for the entity
      const NSMesh::NSFields::EntityComInfo& entity_com_info = it_e.second;
      
      // Check if the processor owns the entity or not
      const int& owning_proc = entity_com_info.GetOwningProcID();
      bool owned = owning_proc == proc_rank;
      
      // Set the number of data to be communicated from this processor to the others
      // This number is by default 2 because we communicate the old and new gid
      unsigned e_size = 2;
      
      // If the processor is the owner of this entity
      // Get the sharing processors' rank first
      std::vector<int> all_non_owning_proc_id = entity_com_info.GetSharingProcID();
      
      // Add the processor's rank for which this entity is universal
      const std::vector<int>& universal_proc_id = entity_com_info.GetUniversalProcID();    
      for( size_t u=0; u<universal_proc_id.size(); ++u ) {
	all_non_owning_proc_id.push_back( universal_proc_id[u] );
      }
      
      if( all_non_owning_proc_id.size() > 0 ) {
	// If the processor is the owner, then it should communicate the new gid to the other processes since it is responsible for setting it
	if( owned ) {
	  // Loop over all the non owning processors
	  for( size_t sp=0; sp<all_non_owning_proc_id.size(); ++sp ) {
	    // Adjust the number of data to be sent to the processor
	    send_size[ all_non_owning_proc_id[sp] ] += 2;
	    
	    // Push back the old gid
	    to_send_gid[ all_non_owning_proc_id[sp] ].push_back( it_e.first );

	    // Push back the new gid 
	    to_send_gid[ all_non_owning_proc_id[sp] ].push_back( p_oldToNewGIDMap.find( it_e.first )->second );
	  }
	}
	else {
	  // Adjust the number of data to be received from the owning processor
	  recv_size[ owning_proc ] += e_size;
	}
      }
    }
    
    // Commnunicate the data between the processors
    std::vector<std::vector<unsigned> > to_recv_gid( com_size, std::vector<unsigned>{} );
    _executeProcComRequestsForGIDs( recv_size, send_size, to_send_gid, to_recv_gid );
    
    // Adjust the gids of the communicated entity data
    for( size_t p=0; p<to_recv_gid.size(); ++p ) {
      const std::vector<unsigned>& gids_vec = to_recv_gid[p];
      
      for( size_t e=0; e<gids_vec.size(); e+=2 ) {
	// Get the entity using the the old gid
	NSMesh::NSEntity::Entity& entity = _getEntity( p_topology, gids_vec[e] );
	
	// Set the new gid of the entity
	entity.SetGid( gids_vec[e+1] );
	
	// Store the old and new gids if the topology of the entity is Topology::tp_VERTEX
	// These are needed later to create a finite element field for the coordinates of the mesh points
	if( p_topology == NSMesh::Topology::tp_VERTEX ) {
	  _newToOldVerticesGIDMap.insert( std::make_pair( gids_vec[e + 1], gids_vec[e] ) );
	}
      }
  }
    
    // Adjust the entity com info map
    // First copy the existing map to a temporary map
    auto copy_entity_com_info_map = _topoToEntityGidProcMap.find( p_topology )->second;

    // Clear the existing map
    _topoToEntityGidProcMap.find( p_topology )->second.clear();
  
    // Loop over the communicated data and add them to the map
    auto& new_entity_com_info_map = _topoToEntityGidProcMap.find( p_topology )->second;
    for( size_t p=0; p<to_recv_gid.size(); ++p ) {
      const std::vector<unsigned>& gids_vec = to_recv_gid[p];
      
      for( size_t e=0; e<gids_vec.size(); e+=2 ) {
	// Add the new entity gid along with its corresponding pointer to the entity com info object
	new_entity_com_info_map.insert( std::make_pair( gids_vec[e+1], std::move( copy_entity_com_info_map.find( gids_vec[e] )->second ) ) );
	
	// Removed the added gid from the copy_entity_com_info_map
	copy_entity_com_info_map.erase( gids_vec[e] );
      }
    }
    
    // Add the remaining entities to the new_entity_com_info_map
    for( auto& it : copy_entity_com_info_map ) {
      new_entity_com_info_map.insert( std::make_pair( p_oldToNewGIDMap.find( it.first )->second, it.second ) );
    }
  }
  
  // Rebuild the gid to lid map of the EntityContainer of all the entities of a specific topology defined on the whole mesh
  _topologyToEntityContainerMap.find( p_topology )->second->RebuildGIDToLIDMap();

  // Rebuild the gid to lid map of the EntityContainer of all the entities of a specific topology defined on each InteriorMeshPart
  for( auto it = _nameToInteriorMeshPartMap.begin(); it != _nameToInteriorMeshPartMap.end(); ++it ) {
    it->second->RebuildGIDToLIDMapPerTopology( p_topology );
  }

  // Rebuild the gid to lid map of the EntityContainer of all the entities of a specific topology defined on each BoundaryMeshPart
  for( auto it = _nameToBoundaryMeshPartMap.begin(); it != _nameToBoundaryMeshPartMap.end(); ++it ) {
    it->second->RebuildGIDToLIDMapPerTopology( p_topology );
  }
}


void
STKMesh::_executeProcComRequestsForGIDs( const std::vector<unsigned>& p_recvSize, const std::vector<unsigned>& p_sendSize, const std::vector<std::vector<unsigned> >& p_toSndGID, std::vector<std::vector<unsigned> >& p_toRecvGID ) {
  // Get the processor rank
  int procRank = _com->GetRank();

  // Get the number of receive requests to be made
  int rcv_cnt = 0;
  for( int rr=0; rr<static_cast<int>(p_recvSize.size()); ++rr ) {
    if( p_recvSize[rr] > 0 ) {
      p_toRecvGID[rr].resize( p_recvSize[rr] );
      rcv_cnt++;
    }
  }

  // Get the number of send requests to be made  
  int snd_cnt = 0;
  for( int sr=0; sr<static_cast<int>(p_sendSize.size()); ++sr ) {
    if( p_sendSize[sr] > 0 ) {
      snd_cnt++;
    }
  }

  // Allocate memory for the send and receive requests
  MPI_Request* rcv_requests = _com->GenerateMPIRequest( rcv_cnt );
  MPI_Status* rcvStatus = _com->GenerateMPIStatus( rcv_cnt );
  MPI_Request* snd_requests = _com->GenerateMPIRequest( snd_cnt );
  MPI_Status* sndStatus = _com->GenerateMPIStatus( snd_cnt );
  
  // Enqueue the recv calls for the field data
  rcv_cnt = 0;
  for( int rr=0; rr<static_cast<int>(p_recvSize.size()); ++rr ) {    
    if( rr != procRank && p_recvSize[rr] > 0 ) {
      _com->IRecv( &p_toRecvGID[rr][0], p_recvSize[ rr ], MPI_UNSIGNED, rr, 2, &rcv_requests[ rcv_cnt ] );
      rcv_cnt++;
    }
  }
  
  // Enqueue the send calls for the field data
  snd_cnt = 0;
  for( int sr=0; sr<static_cast<int>(p_sendSize.size()); ++sr ) {    
    if( sr != procRank && p_sendSize[sr] > 0 ) {
      _com->ISend( &p_toSndGID[sr][0], p_sendSize[sr], MPI_UNSIGNED, sr, 2, &snd_requests[ snd_cnt ] );
      snd_cnt++;
    }
  }

  // Wait for the send calls to be completed
  if( snd_cnt > 0 ) {
    MPI_Waitall( snd_cnt, snd_requests, sndStatus );
  }

  // Wait for the receive calls to be completed
  if( rcv_cnt > 0 ) {
    MPI_Waitall( rcv_cnt, rcv_requests, rcvStatus );
  }

  // Clear the memory allocated for the send and receive requests
  rcv_requests = _com->ClearMPIRequest( rcv_requests );
  rcvStatus = _com->ClearMPIStatus( rcvStatus );
  snd_requests = _com->ClearMPIRequest( snd_requests );
  sndStatus = _com->ClearMPIStatus( sndStatus );
}

NSMesh::ElementTopology
STKMesh::_getEltypeFromCellTopology( const shards::CellTopology& p_cellTopology ) const
{
  if ( p_cellTopology == shards::getCellTopologyData<shards::Line<2> >() ) {
    return NSMesh::ElementTopology::et_LINE2;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Line<3> >() ) {
    return NSMesh::ElementTopology::et_LINE3;
  }
  if ( p_cellTopology == shards::getCellTopologyData<shards::Triangle<3> >() ) {
    return NSMesh::ElementTopology::et_TRI3;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Triangle<6> >() ) {
    return NSMesh::ElementTopology::et_TRI6;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Quadrilateral<4> >() ) {
    return NSMesh::ElementTopology::et_QUAD4;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Quadrilateral<8> >() ) {
    return NSMesh::ElementTopology::et_QUAD8;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Tetrahedron<4> >() ) {
    return NSMesh::ElementTopology::et_TET4;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Tetrahedron<10> >() ) {
    return NSMesh::ElementTopology::et_TET10;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Hexahedron<8> >() ) {
    return NSMesh::ElementTopology::et_HEX8;
  }
  else if ( p_cellTopology == shards::getCellTopologyData<shards::Hexahedron<20> >() ) {
    return NSMesh::ElementTopology::et_HEX20;
  }
  else {
    NSExceptions::InvalidArgument( true, "The mesh element topology written in the mesh file is not supported by the code" );
    return NSMesh::ElementTopology::et_TRI3;
  }
}
    
} // namespace STKMesh
} // namespace Mesh
