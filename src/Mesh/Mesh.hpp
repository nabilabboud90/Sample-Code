#ifndef MESH_HPP_
#define MESH_HPP_

/*
 * @brief Abstract class to represent a mesh.
 * For now this class reads the mesh using the stk package in trilinos
 * It then extracts mesh information such as connectivity and gids from the stk mesh data strcutures and stores it in its own data structure
 * This class also holds a pointer a mesh query class that allows the storing of a variation of information about the mesh
 * In particular depending on the mesh query class chose, the code will store explicitly some connectivity information and will rely on implicit routines to obtain other connectivity information
 * For more details about the different information that can be stored about the mesh the reader is refered to
 * "Mesh data structure selection for mesh generation and FEA application", Rao Garimella
 */

#include<string>
#include<memory>
#include<map>
#include<vector>

#include "Utility/Constants.hpp"
#include "Mesh/MeshEnum.hpp"
#include "Entity/EntityContainer.hpp"
#include "Entity/EntityStats.hpp"
#include "Entity/Entity.hpp"
#include "MeshPart/InteriorMeshPart.hpp"
#include "MeshPart/BoundaryMeshPart.hpp"
#include "MeshUtility/EntityComInfo.hpp"

#include <Kokkos_Core.hpp>

namespace NSFactory { template<typename Product, typename Identifier, typename... Argument> class Factory; }
namespace NSParser { namespace NSInput { class Datafile; } }
namespace NSCommunicator { namespace NSMPIComWrapper { class MPIComWrapper; } }

namespace NSMesh {
  
class MeshQuery;

namespace NSFields { class BasicMeshField; }

namespace NSEntity { template<typename T> class Bucket; }

class Mesh {
public:
  typedef NSFactory::Factory<Mesh*, std::string, void> factory_mesh_Type; /* *< Factory to chose among the different mesh input files supported in the code (stk, gmesh, etc) */

  /*
   * Default constructor
   */
  explicit Mesh();

  /*
   * Copy constructor not implemented
   @param mesh - Mesh object used to copy from
   */
  Mesh( const Mesh& p_mesh ) = delete;

  /*
   * Equal operator not implemented
   @param mesh - Mesh object used to copy from
  */
  Mesh& operator=( const Mesh& p_mesh ) = delete;

  /*
   * Setup function that parses all the parameters needed by this class from the input file
   @param datafile - reference to the input file
   @param sectionPath - a string representing the path to the section in the input file where the mesh parameters are specified
   @param sectionName - name of the section in the input file that holds all the information about the mesh parameters
  */
  virtual void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Initialization function
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  virtual void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com );
  
  /*
   * Function to initialize the coordinates field
   @param bucket - Bucket of elements
   @param coordinates - Vector that will hold the values of the coordinates of the vertices of the elements
  */
  virtual void PopulateCoordinatesField( const NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*>& p_bucket, const Kokkos::View<double***>& p_coordinates ) const = 0;

  /*
   * @brief The way the data stored is as follows
   * There are two different type of data that are stored by the mesh which are global ids of entities (vertex, edges, faces, and elements) and connectivity between those entities
   * The global ids are stored as unsigned int and the connectivities are stored as pointers
   * Each entity in the mesh is represented by an Entity class that holds an unsigned int representing the global ids of the entity and a map of pointers to the connected entities of different topology such as element, faces, etc
   * The amount of data stored in the map about the connectivity is based on the mesh query type that is specified in the input file
   * If no query type is specified a default one is used, check the code to determine the default type as it may change based on performance analysis
   * For each entity topology (i.e, vertex, element, face, edge) there exists an entity container that will hold all the entities of a given topology
   * The entity container class serves as a classifier of the different entities of the mesh based on their topology
   * Within each entity container, there exists a vector of entities that holds the all the entities of a given topology (for instance all the vertices)
   * The entities in this vector are stored in a particular order as follow
   * First all the owned entities of a given topology are stored (for instance all the owned vertices)
   * Second all the shared entities of a given topology are stored (for instance all the shared vertices)
   * Third all the universsal entities of a given topology are stored (for instance all the universal vertex)
   * To keep track of the indexes of the entity vector where entities of a given ownership type ends and entities of a different ownership type starts an additional information is stored in the entity container
   * This additional information is stored as a map from an ownership type to a vector of pairs of long unsigned int where those pair of long unsigned int represent the start and end index respectively of those entities in the entity vector
   * The reason i use a vector of pair of long unsigned and not just a pair of long unsigned int is because if we later decide to add entities to the mesh (for instance mesh refinement) then we would like to add those entities to the end of the entities vector in order to avoid memory reallocation as much as possible, hence, in this case the entities of a given ownership type won't be contiguous anymore in the entity vector
   */

  /*
   * Function that returns pointers to entities given the gids of those entities
   @param buckets - Vector of Bucket where each bucket represents a container of the entities to be returned
   @param entityGids - A vector holding the gids of the entities whose pointers will be returned
   @param topology - The topology of the entities whose gid is given in the entityGids
   @param capacity - Maximum number of entities that can be supported by a bucket. This number is set by default but can be changed if needed
  */
  void GetBucketsByGid( std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> >& p_buckets, const std::vector<unsigned>& p_entityGids,
			const NSMesh::Topology& p_topology, unsigned p_capacity = BUCKETSIZE ) const;
  
  /*
   * Get function that returns the entities of a certain type of a given ownership by the processor (owned, shared, universal)
   * Note that a bucket will hold a pointer a given entity in memory and a size from which the number of entites that will be represented by this bucket can be determined
   @param buckets - Vector of Bucket where each bucket represents a container of the entities to be returned
   @param ownership - Is a vector of the type of entity ownership we are asking for (owned, shared, universal or a union of them)
   @param topology - Entity topology we are asking for (vertex, element, edge, face)
   @param capacity - Maximum number of entities that can be supported by a bucket. This number is set by default but can be changed if needed
   */
  void GetBuckets( std::vector<NSMesh::NSEntity::Bucket<std::unique_ptr<NSMesh::NSEntity::Entity> > >& p_buckets,
		   const std::vector<NSMesh::Ownership>& p_ownership, const NSMesh::Topology& p_topology, unsigned p_capacity = BUCKETSIZE ) const;

  /*
   * Get function that returns the global ids (gids) of a certain entity type of a given ownership by the processor (owned, shared, universal)
   * Note that a bucket will hold a pointer a given entity in memory and a size from which the number of entites that will be represented by this bucket can be determined
   @param buckets - Vector of Bucket where each bucket represents a container of the entities to be returned
   @param ownership - Is a vector of the type of entity ownership we are asking for (owned, shared, universal or a union of them)
   @param topology - Entity topology we are asking for (vertex, element, edge, face)
   @param interiorPartName - Name of the interior mesh part from which to get the entities
   @param capacity - Maximum number of entities that can be supported by a bucket. This number is set by default but can be changed if needed
   */
  void GetBucketsByInteriorPart( std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> >& p_buckets, const std::vector<NSMesh::Ownership>& p_ownership, const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, unsigned p_capacity = BUCKETSIZE ) const;

  /*
   * Get function that returns the global ids (gids) of a certain entity type of a given ownership by the processor (owned, shared, universal)
   * Note that a bucket will hold a pointer a given entity in memory and a size from which the number of entites that will be represented by this bucket can be determined
   @param buckets - Vector of Bucket where each bucket represents a container of the entities to be returned
   @param ownership - Is a vector of the type of entity ownership we are asking for (owned, shared, universal or a union of them)
   @param topology - Entity topology we are asking for (vertex, element, edge, face)
   @param boundaryPartName - Name of the boundary mesh part from which to get the entities
   @param interiorPartName - Name of the interior mesh part that is attached to the boundary mesh part from whom we are gathering the entities
   @param capacity - Maximum number of entities that can be supported by a bucket. This number is set by default but can be changed if needed
   */
  void GetBucketsByBoundaryPart( std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> >& p_buckets, const std::vector<NSMesh::Ownership>& p_ownership, const NSMesh::Topology& p_topology, const std::string& p_boundaryPartName, const std::string& p_interiorPartName, unsigned p_capacity = BUCKETSIZE ) const;
  
  /*
   * Get function that returns an entity based on its topology and gid
   @param topology - Parameter specifying the topology of the entity we are asking for (vertex, element, edge, face)
   @param gid - Global id of the entity to be returned
   @output Entity that has topology and gid
  */
  const NSMesh::NSEntity::Entity& GetEntity( const NSMesh::Topology& p_topology, const unsigned& p_gid ) const;

  /*
   * Function that returns the topology of an entity
   @param gid - The global id of the entity
   @output The topology of the entity whose gid is given
  */
  NSMesh::Topology GetEntityTopology( const unsigned& p_gid ) const;
  
  /*
   * This is a connectivity related routine. It returns a vector of pointers to entities of a given topology that are connected to entity
   @param topology - Topology of the connected entities we are asking for
   @param entity - Entity that we want to query its connected entities with the specified topology
   @output connectedEntitiesVector - Vector of pointers to the entities of topology and connected to entity
  */
  void Begin( const NSMesh::Topology& p_toTopology, const NSMesh::NSEntity::Entity& p_entity, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVector ) const;

  /*
   * Getter function that returns the name of the mesh
   @output Name of the mesh
  */
  const std::string& GetName() const;
  
  /*
   * Getter function of the number of an entity for a specified set of ownership type per processor
   @param ownership - Vector of ownership types for which to return the total number of a given entity
   @param topology - Topology of the entity we are querying
   @output Total number of a given entity of the specified onwership type per processor
  */
  unsigned GetNumberOfEntity( const NSMesh::Topology& p_topology, const std::vector<NSMesh::Ownership>& p_ownership ) const;

  /*
   * Getter function of the number of an entity for a specified set of ownership type per processor for a set of parts
   @param ownership - Vector of ownership types for which to return the total number of a given entity
   @param topology - Topology of the entity we are querying
   @param - List of parts on which the number of entities is needed
   @output Total number of a given entity of the specified onwership type per processor
  */
  unsigned GetNumberOfEntityPerParts( const NSMesh::Topology& p_topology, const std::vector<NSMesh::Ownership>& p_ownership, const std::vector<std::string>& p_parts ) const;
  
  /*
   * Getter function of the total number of an entity on the whole mesh
   @param topology - Topology of the entity we are querying
   @output Total number of a given entity on the whole mesh
  */
  const size_t& GetGlobalNumberOfEntity( const NSMesh::Topology& p_topology ) const;

  /*
   * Function that deletes the third party library mesh data structure if used
   */
  virtual void ClearTPLMesh() = 0;

  /*
   * Getter function of spatial dimension of the problem
   @output Spatial dimension of the problem
   */
  int SpatialDimension() const;

  /*
   * Getter function of the communicator that is used to communciate data defined on the mesh
   @output The wrapper class around the mpi communicator
  */
  const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& GetMPIComWrapper() const;

  /*
   * Getter function of the data structure holding the entities at the interface between processors
   @output Data structure holding the data about the entities at the interface between processors
  */
  const std::map<NSMesh::Topology, std::map<unsigned, NSMesh::NSFields::EntityComInfo> >& GetProcInterfaceEntityInfo() const;

  /*
   * Function that checks whether a mesh part with a specified name is defined
   @param name - Name of the mesh part to check for
   @output Boolean specifying whether the mesh part has been defined
  */
  bool IsMeshPartDefined( const std::string& p_name ) const;

  /*
   * Getter function that returns the topology of the closure of the domain/element of the mesh
   @output Closure topology
  */
  NSMesh::Topology GetClosureTopology() const;

  /*
   * Function that checks whether an entity is on the interface between processors
   * For an entity to be considered on the processors interface it must be either shared or universal with respect to another processor
   @param gid - Global id of the entity being checked
   @output Boolean specifying whether the entity is on the interface between processes or not
  */
  bool IsOnProcessorInterface( unsigned p_gid ) const;

  /*
   * Function that checks whether an entity is owned by a processor
   * To be owned, the entity must either not be on the interface between processors or if it is, the id of the calling processor must match the owner rank assigned to the entity by the mesh partitioner
   @param gid - The global id of the entity being checked
   @output A boolean specifying if the entity is owned or not
   */
  bool IsOwned( unsigned p_gid ) const;

  /*
   * Function that returns the id of the owning processor of the entity
   @param gid - The global id of the entity whose owning processor will be returned
   @output The id of the owning processor
  */
  unsigned GetOwningProcessor( unsigned p_gid ) const;
  
  /*
   * Getter function of the map from an entity gids to a pointer to the object holding the data about the owning, sharing, and universal processes
   * This list includes shared entities only and not universal
   @param topology - Topology of the interface entities whose gids will be returned
  */
  const std::map<unsigned, NSMesh::NSFields::EntityComInfo>& GetProcessorInterfaceEntities( const NSMesh::Topology& p_topology ) const;

  /*
   * Getter function of the name of all the mesh part of a given topology
   @param topology - Topology of the primary entities of the mesh parts requested
   @output A vector holding the name of all the mesh parts of a given topology defined on this mesh
  */
  std::vector<std::string> GetAllMeshPartsName( const NSMesh::Topology& p_topology ) const;

  /*
   * Getter function of the mesh element topology for a given mesh part
   @param partName - Name of the mesh part whose element topology will be returned
   @output The topology of the element that define the mesh part
  */
  NSMesh::ElementTopology GetPartElementTopology( const std::string& p_partName ) const;

  /*
   * Getter function of the number of mesh verticess per element for a given mesh part
   * This function assumes that all the mesh elements are of the same topology
   @param partName - The name of the mesh part whose number of mesh vertices per element will be returned
   @output Number of mesh vertices per element of the mesh
  */
  int GetNumberOfVerticesPerElement( const std::string& p_partName ) const;

  /*
   * Getter function of a boundary mesh part
   @param partName - The name of the boundary mesh part to be returned
   @output A pointer to the mesh part
  */
  const NSMesh::BoundaryMeshPart* GetBoundaryMeshPart( const std::string& p_partName ) const;

  /*
   * Getter function of an interior mesh part
   @param partName - The name of the interior mesh part to be returned
   @output A pointer to the mesh part
  */
  const NSMesh::InteriorMeshPart* GetInteriorMeshPart( const std::string& p_partName ) const;

  /*
   * Function that sorts a set of entities based on their topology
   @param entitiesGid - A vector holding the global ids of all the entities to be sorted
   @param sortedEntitiesGid - A 2D vector holding the sorted entities gids
  */
  void SortByTopology( const std::vector<unsigned>& p_entitiesGid, std::vector<std::vector<unsigned> >& p_sortedEntitiesGid ) const;
  
private:
  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - pointer to the input file
   @param section_path - a string representing the path to the section in the input file where the mesh parameters are specified
   @param section_name - name of the section in the input file that holds all the information about the mesh parameters
  */
  void _setupMeshQuery( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - pointer to the input file
   @param section_path - a string representing the path to the section in the input file where the mesh parameters are specified
   @param section_name - name of the section in the input file that holds all the information about the mesh parameters
  */
  void _setupMeshParts( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Function that initializes the _topology_to_entiy_container_map
   * The call to this function will not fill the entity containers, it will just create and set them in memory
   */
  void _initializeTopologyToEntityContainerMap();

  /*
   * Function that initialize the mesh parts
   */
  virtual void _initializeMeshParts() = 0;
  
  /*
   * Function to populate the entity containers stored in the _topologyToEntityContainerMap with the corresponding entities
   * The implementation of this function depends on the underlying mesh file being read (different if using stk or using gmsh for instance)
   */
  virtual void _populateTopologyToEntityContainerMap() = 0;

  /*
   * Function to populate the connectivity relations of the entities stored in the entity containers of the _topologyToEntityContainerMap 
   * The implementation of this function depends on the underlying mesh file being read (different if using stk or using gmsh for instance)
   */
  virtual void _populateEntityConnectivityInformation() = 0;

  /*
   * Function to populate the mesh parts with there respective entities
   * The mesh part will hold the gid of those entities and not the entity object itself for memory efficiency purposes since those entities are duplicated and already stored somewhere else
   */
  virtual void _populateMeshPartsEntities() = 0;

  /*
   * Function that populates the communication information about the owned, shared and universal entities
   */
  virtual void _populateEntityComInfo() = 0;

  /*
   * Function to print a summary of the number of entities defining the mesh
   */
  void _printMeshSummary() const;

  /*
   * Function that populates the buckets with entities of a certain topology
   @param buckets - Vector of Bucket where each bucket represents a container of the entities to be returned
   @param indexesVec - Vector holding the indices of entities in the vector where they are stored
   @param entityVec - Vector holding pointers to entities of a specified topology and ownership
   @param capacity - Maximum number of entities that can be supported by a bucket. This number is set by default but can be changed if needed
  */
  template<typename T1, typename T2>
  void _populateBuckets( std::vector<NSMesh::NSEntity::Bucket<T1> >& p_buckets,
			 const std::vector<std::pair<long unsigned, long unsigned> >& p_indexesVec,
			 const std::vector<T2>& p_entityVec, unsigned p_capacity ) const;
  
protected:
  /*
   * Get function that returns an entity based on its topology and gid
   @param topology - Parameter specifying the topology of the entity we are asking for (vertex, element, edge, face)
   @param gid - Global id of the entity to be returned
   @output Entity that has topology and gid
  */
  NSMesh::NSEntity::Entity& _getEntity( const NSMesh::Topology& p_topology, const unsigned& p_gid );

  /*
   * Renumber entities
   */
  virtual void _renumberEntities() = 0;
  
protected:
  std::string _name = {}; /**<Mesh name*/
  std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper> _com = {}; /**<Communicator responsible for doing any communication on the mesh*/
  int _ndim; /**<Number of dimensions of the problem we are solving*/
  std::map<NSMesh::Topology, std::unique_ptr<NSMesh::NSEntity::EntityContainer<std::unique_ptr<NSMesh::NSEntity::Entity> > > > _topologyToEntityContainerMap = {}; /**<Map that returns given a topology, all the entities in the mesh with that topology. Those entities are stored in an entity container*/
  std::unique_ptr<NSMesh::MeshQuery> _meshQuery; /**<Pointer to the mesh query class that holds all the protocols regarding the connectivity information to be explicitly and implicitly stored*/
  std::map<std::string, std::unique_ptr<NSMesh::InteriorMeshPart> > _nameToInteriorMeshPartMap = {}; /**<Map of the interior mesh part name to the pointer to the interior mesh part*/
  std::map<std::string, std::unique_ptr<NSMesh::BoundaryMeshPart> > _nameToBoundaryMeshPartMap = {}; /**<Map of the boundary mesh part name to the pointer to the boundary mesh part*/
  std::map<NSMesh::Topology, std::map<unsigned, NSMesh::NSFields::EntityComInfo> > _topoToEntityGidProcMap = {}; /**<This data structure holds data about entities defined at the interface between processors. The data consists of the processor that owns, shares or have access to this entity*/
  bool _renumberMeshEntities = true; /**<Boolean to specify whether to renumber the entities after they have been read from the mesh file. For details about the new renumbering scheme check the method RenumberEntities */
  bool _createAllSideEntities = {}; /**<Boolean specifying whether to create edges and faces for all the elments of the domain or only for for those elements on the boundary/interface of the domain*/
  std::map<unsigned, unsigned> _newToOldVerticesGIDMap = {}; /**<This map stores the new gids and the old gids of the vertices if the vertices are renumbered*/
};
  
} // namespace NSMesh

#endif
