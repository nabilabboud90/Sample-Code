#ifndef STKMESH_HPP_
#define STKMESH_HPP_

/*
 * @brief STK mesh wrapper that reads and extracts all the needed mesh information using the stk library
 * Once the information are read, they will be stored in data strucutres that are members of the parent class mesh
 * Once all the information is stored, the stk data structures are destroyed to free up the memory
 */

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <stk_mesh/base/Types.hpp>
#include <Teuchos_RCP.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#pragma GCC diagnostic pop

#include "Mesh/Mesh.hpp"
#include "Factory/Factory.hpp"

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }

namespace NSMesh { namespace NSFields { class BasicMeshField; } }

namespace NSMesh { namespace NSSTKMesh {

class STKMesh : public Mesh {
public:
  /*
   * Default constructor
   */
  explicit STKMesh();

  /*
   * Copy constructor not implemented
   @param mesh - STK mesh to copy from
  */
  STKMesh( const STKMesh& p_mesh ) = delete;

  /*
   * Equal operator not implemented
   @param mesh - STK mesh to set equal to
  */
  STKMesh& operator=( const STKMesh& p_mesh ) = delete;

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - pointer to the input file
   @param section_path - a string representing the path to the section in the input file where the stk mesh parameters are specified
   @param section_name - name of the section in the input file that holds all the information about the stk mesh parameters
  */  
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) override;

  /*
   * Initialization function
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) override;

  /*
   * Function that deletes the third party library mesh data structure if used
   */
  void ClearTPLMesh();

  /*
   * Function to initialize the coordinates field
   @param bucket - Bucket of elements
   @param coordinates - Vector that will hold the values of the coordinates of the vertices of the elements
  */
  void PopulateCoordinatesField( const NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*>& p_bucket, const Kokkos::View<double***>& p_coordinates ) const override;

private:
  stk::mesh::MetaData* _metaData = {}; /**<Class in stk that holds information about the mesh*/
  std::unique_ptr<stk::mesh::BulkData> _bulkData = {}; /**<Class in stk that holds information about the mesh*/
  std::unique_ptr<stk::io::StkMeshIoBroker> _meshReader = {}; /**<Class in stk that holds information about the mesh*/
  std::map<NSMesh::Topology, stk::mesh::EntityRank> _topologyToSTKMeshTopology = {}; /**<Map from the topology enumerations in this code to the topology enumerations in stk*/
  std::map<std::string, stk::mesh::Part*> _partNameToStkMeshPart = {}; /**<Map from the stk topology to a vector of pointers to stk mesh parts with that topology*/
  std::map<std::string, std::pair<NSMesh::ElementTopology, int> > _partNameToMeshCellTopologyInfo = {}; /**<Map from the part name to information about the topology of the cells defining this part. The assumption is that all cells within one part have the same topology*/
  
private:
  /*
   * Function that initialize the mesh parts
   */
  void _initializeMeshParts() override;
  
  /*
   * Function to populate the entity containers stored in the _topologyToEntityContainerMap with the corresponding entities
   */
  void _populateTopologyToEntityContainerMap() override;

  /*
   * Function to populate the entity container by creating the entities given a topology and ownership type
   @param rank - STK topology enum that specifies the topology of the entities to be created to populate the entity container with
   @param selector - STK ownership type to specify the entity ownership type to be created to populate the entity container with
   @param topology - Topology enum that specifies the topology of the entities to be created to populate the entity container with
   @param offset - Starting index in the entities vector of the entity of a given ownerhips type
   */
  void _populateEntityContainer( const stk::mesh::EntityRank& p_rank, const stk::mesh::Selector& p_selector, const NSMesh::Topology& p_topology, unsigned p_offset = 0 );

  /*
   * Function that returns the total number of entities in the mesh of a given topology and ownership type
   @param rank - STK topology enum
   @param selector - STK ownership type
   @output Total number of entities in the mesh of a given topology and ownership type
  */
  size_t _getTotalNumEntities( const stk::mesh::EntityRank& p_rank, const stk::mesh::Selector& p_selector ) const;

  /*
   * Function to populate an entity container of vertices
   @param rank - STK topology enum that specifies the topology of the entities to be created to populate the entity container with
   @param topology - Topology enum that specifies the topology of the entities to be created to populate the entity container with
   */
  void _populateEntityInfo( const stk::mesh::EntityRank& p_rank, const NSMesh::Topology& p_topology );

  /*
   * Function to populate the connectivity relations of the entities stored in the entity containers of the _topologyToEntityContainerMap 
   */
  void _populateEntityConnectivityInformation() override;

  /*
   * Function to populate the connectivity relations between two entities given the topology of the source entity, the topology of the connected entity, and the ownership type
   @param from_entity_topology - Topology enum of the entity for which we would like to store the connectivity information
   @param to_entity_topology - Topology enum of the connected entities to the entity we would like to store the information for
   @param selector - STK ownership enum of the entities we would like to store the connectivity information for
  */
  void _populateEntityToEntityConnectivity( const NSMesh::Topology& p_fromEntityTopology, const NSMesh::Topology& p_toEntityTopology, const stk::mesh::Selector& p_selector );

  /*
   * Function that creates the stk mesh parts that will be used to define the mesh parts
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
   */
  void _createSTKMeshParts( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com );

  /*
   * Function to populate the mesh parts with there respective entities
   * The mesh part will hold the gid of those entities and not the entity object itself for memory efficiency purposes since those entities are duplicated and already stored somewhere else
   */
  void _populateMeshPartsEntities() override;
  
  /*
   * Function to populate the interior mesh parts with there respective entities
   * The mesh part will hold the gid of those entities and not the entity object itself for memory efficiency purposes since those entities are duplicated and already stored somewhere else
   */
  void _populateInteriorMeshPartsEntities();

  /*
   * Function that populates the entities of a specific topology of an interior mesh part
   @param meshPart - A pointer to the interior mesh part to be populated
   @param stkMeshPart - A pointer to the interior mesh part of type stk
   @param topology - The topology of the entities belonging to the mesh part whose information will be stored in the mesh part
  */
  void _populateInteriorMeshPartsEntitiesPerTopology( NSMesh::InteriorMeshPart* const p_meshPart, const stk::mesh::Part& p_stkMeshPart,
						      const NSMesh::Topology& p_topology );
  
  /*
   * Function to populate the vector of pointers to entities of a given mesh part
   @param name - Name of the mesh part for which to populate the vector of pointers
   @param stk_topology - STK topology of the entities for which to populate the vector of pointers
   @param selector - Subset of the mesh from which to collect the entities
   @param topology - Mesh topology of the entities for which to populate the vector of pointers
   @param offset - Starting position in the entity pointers vector of a given mesh part from where to start to populate the entity pointer
  */
  void _populateInteriorMeshPartEntityVector( const std::string& p_name, const stk::mesh::EntityRank& p_stkTopology, const stk::mesh::Selector& p_selector, const NSMesh::Topology& p_topology, unsigned p_offset );

  /*
   * Function that populates the boundary mesh parts
   */
  void _populateBoundaryMeshPartsEntities();

  /*
   * Function that populates a given mesh part with information about the closure entities and the elements attached to those closure entities
   @param meshPart - A pointer to the boundary mesh part
   @param stkMeshPart - A pointer to the boundary mesh part of type stk
  */
  void _populateBoundaryMeshPartsEntitiesForAllTopologies( NSMesh::BoundaryMeshPart* const p_meshPart, const stk::mesh::Part& p_stkMeshPart );

  /*
   * A function that seggregates the clsoure entities and connected elements to those entities based on the ids of the domains to which those closure entities are connected
   @param stkTopology - Stk topology of the closure entities
   @param selector - A selector referencing the entities being processed
  */
  void _populateBoundaryMeshPartEntityVector( const stk::mesh::EntityRank& p_stkTopology, const stk::mesh::Selector& p_selector,
					      std::map<std::string, std::vector<unsigned> >& p_interiorPartNameToClosureEntityGid,
					      std::map<std::string, std::vector<unsigned> >& p_interiorPartNameToElementGid );

  /*
   * Function that stores the data about entities of a boundary mesh part of a certain topology that are attached to given interior mesh part
   @param interiorPartName - The name of the interior mesh part to whom the boundary mesh part is attached
   @param boundaryMeshPart - A reference to the boundary mesh part object
   @param ownedEntityGid - A vector holding the gids of the owned entities
   @param sharedEntityGid - A vector holding the gids of the shared entities
   @param universalEntityGid - A vector holding the gids of the universal entities
   @param topology - The topology of the entities being stored
  */
  void _completeBoundaryMeshPartPerTopology( const std::string& p_interiorPartName, NSMesh::BoundaryMeshPart& p_boundaryMeshPart,
					     const std::vector<unsigned>& p_ownedEntityGid, const std::vector<unsigned>& p_sharedEntityGid,
					     const std::vector<unsigned>& p_universalEntityGid, const NSMesh::Topology& p_topology );

  /*
   * Add entities of a certain topology to the boundary mesh part
   @param entityGids - Vector holding the gids of the entities to be added
   @param meshPart - A reference to the boundary mesh part
   @param interiorPartName - The name of the interior mesh part to whom the boundary mesh part is attached
   @param offset - Offset used to set the local ids of the entities added to the boundary mesh part
   @param topology - Topology of the entities being added
  */
  void _addEntityGidToBoundaryMeshPart( const std::vector<unsigned>& p_entityGids, NSMesh::BoundaryMeshPart& p_meshPart,
					const std::string& p_interiorPartName, unsigned p_offset, const NSMesh::Topology& p_topology );

  /*
   * Function that synchronizes the size of temporary maps holding data about entities of a specific topology on a processor
   @param ownedData - Map holding the owned data of a specific topology
   @param sharedData - Map holding the shared data of a specific topology
   @param universalData - Map holding the universal data of a specific topology
  */
  void _synchronizeBoundaryMeshMapsPerTopologyPerProcessor( std::map<std::string, std::vector<unsigned> >& p_ownedData,
							    std::map<std::string, std::vector<unsigned> >& p_sharedData,
							    std::map<std::string, std::vector<unsigned> >& p_universalData ) const;

  /*
   * Function that synchronizes the size of temporary maps holding data about entities of all topologies and across all processors
   @param closureEntityOwnedData - Map holding the owned data about closure entities
   @param closureEntitySharedData - Map holding the shared data about closure entities
   @param closureEntityUniversalData - Map holding the universal data about closure entities
   @param closureEntityOwnedData - Map holding the owned data about elements
   @param closureEntitySharedData - Map holding the shared data about elements
   @param closureEntityUniversalData - Map holding the universal data about elements
  */
  void _synchronizeBoundaryMeshMapsPerTopologyPerProcessor( std::map<std::string, std::vector<unsigned> >& p_closureEntityOwnedData,
							    std::map<std::string, std::vector<unsigned> >& p_closureEntitySharedData,
							    std::map<std::string, std::vector<unsigned> >& p_closureEntityUniversalData,
							    std::map<std::string, std::vector<unsigned> >& p_elementOwnedData,
							    std::map<std::string, std::vector<unsigned> >& p_elementSharedData,
							    std::map<std::string, std::vector<unsigned> >& p_elementUniversalData ) const;   
  
  /*
   * Function that populates the communication information about the owned, shared and universal entities
   */
  void _populateEntityComInfo() override;

  /*
   * Function that populates the communication information about the owned, shared and universal entities of a given topology
   @param topology - Entity type for which to collect the info
   @param entities - Vector containing the entite of a given topology for which to collect the info
   */  
  void _populateEntityComInfoForRank( const NSMesh::Topology& p_topology, const std::vector<stk::mesh::Entity>& p_entities );

  /*
   * Renumber entities
   */
  void _renumberEntities() override;
  
  /*
   * Function that renumbers the gids of entities
   @param topology - Topology type whose entities will be renumbered
   @param rank - Rank of the entity that will be renumbered
   @param oldTONewGIDMap - The map that holds the old and new gids of the renumbered entities
   @param offset - A number to offset the starting range of the renumbering of entities on a given processor
  */
  void _renumberEntitiesImpl( const NSMesh::Topology& p_topology, const stk::mesh::EntityRank& p_rank, std::map<unsigned, unsigned>& p_oldToNewGIDMap, size_t p_offset );

  /*
   * Function that communicates the new gids of the shared and universal entities to the other processors
   @param topology - Topology type whose entities will be renumbered
   @param oldToNewGIDMap - The map that holds the old and new gids of the renumbered entities
  */
  void _communicateNewGIDs( const NSMesh::Topology& p_topology, const std::map<unsigned, unsigned>& p_oldToNewGIDMap );

  /*
   * Function to communicate the old and new gids between processors
   @param recv_size - Vector specifying the size of the buffer to be received from the different processors
   @param send_size - Vector specifying the size of the buffer to be sent to the different processors
   @param toSndGID - Buffer of gids to be sent. The assumption is that the gids are held in pairs of old/new
   @param toRecvGID - Buffer that will hold the received gids
  */
  void _executeProcComRequestsForGIDs( const std::vector<unsigned>& p_recvSize, const std::vector<unsigned>& p_sendSize, const std::vector<std::vector<unsigned> >& p_toSndGID, std::vector<std::vector<unsigned> >& p_toRecvGID );

  /*
   * Function that will create the stk edges
   */
  void _createSTKEdges();

  /*
   * Function that will create the stk faces
   */
  void _createSTKFaces();

  /*
   * Function that returns the topology of an element based on the stk cell topology object
   @param cellTopology - Object from the stk library that holds information about the topology of a mesh cell
   @output An enum representing the topology of a mesh element. Note this has nothing to do with the finite element
  */
  NSMesh::ElementTopology _getEltypeFromCellTopology( const shards::CellTopology& p_cellTopology ) const;

  /*
   * Function that communicates the mesh cell topology info about each part of the mesh
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
   */
  void _communicateMeshCellTopologyInfo( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com );

  /*
   * Function that creates the stk mesh parts of a specific type
   @param nameToMeshPartMap - A map from the name of a mesh part to a pointer to the object representing the mesh part of a specific type
  */
  template<typename T>
  void _createSTKMeshParts( const T& p_nameToMeshPartMap );
};

/*
* Create function used by the factory to create this class
@output A mesh whose static type is Mesh and whose dynamic type is STKMesh
*/
Mesh* createSTKMesh();

namespace {
  static bool register_STKMesh = Mesh::factory_mesh_Type::Register( "STKMESH", &createSTKMesh );
}

} // namespace STKMesh
} // namespace Mesh

#endif
