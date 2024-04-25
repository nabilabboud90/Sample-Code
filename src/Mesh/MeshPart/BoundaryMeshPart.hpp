#ifndef BOUNDARYMESHPART_HPP_
#define BOUNDARYMESHPART_HPP_

#include<string>
#include<memory>

#include "Types/Types.hpp"
#include "Mesh/Entity/EntityContainer.hpp"

namespace NSParser { namespace NSInput { class Datafile; } }

namespace NSMesh {
  
class BoundaryMeshPart {
public:
  /*
   * Constructor
   */
  explicit BoundaryMeshPart();

  /*
   * Copy constructor not implemented
   @param mP - Mesh part to copy from
  */
  BoundaryMeshPart( const BoundaryMeshPart& p_mP ) = delete;

  /*
   * Equal operator not implemented
   @param mP - Mesh part to set equal to
  */
  BoundaryMeshPart& operator=( const BoundaryMeshPart& p_mP ) = delete;

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - pointer to the input file
   @param section_path - a string representing the path to the section in the input file where the mesh parameters are specified
   @param section_name - name of the section in the input file that holds all the information about the mesh parameters
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Initialize the boundary mesh part
   @param ndim - The dimension of the problem
  */
  void Initialize( int p_ndim );
  
  /*
   * Getter function of the mesh part name
   @output The name of the mesh part
  */
  const std::string& GetName() const;

  /*
   * Getter function that returns the flag of the mesh part
   @output Flag of the mesh part
  */
  const int& GetId() const;

  /*
   * Getter function that returns the primary topology of the mesh part
   @output Topology of the mesh part
  */
  const NSMesh::Topology& GetPrimaryTopology() const;

  /*
   * Function that sets the information about the topology of the elements that make up the part
   @param numVerticesPerElement - The number of vertices per element
   @param elementTopology - An enum that specifies the shape of the element
  */
  void SetElementTopologyInfo( int p_numVerticesPerElement, const NSMesh::ElementTopology& p_elementTopology );
  
  /*
   * Function to allocate memory to store the entities of a topology defined on this mesh part
   @param topology - Topology of the entities for which to allocate memory
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @param size - Number of entities to allocated memory for
  */
  void AllocateMemoryPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, size_t p_size );

  /*
   * Function that sets statistical data, of a given ownership, about the entities of a specific topology defined on this mesh part
   @param topology - Topology of the entities for which we are setting the statistical data
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @param ownership - Ownership type of the entities for which we are storing the information
   @param data - Number of entities of a given ownership and of a given topology
  */
  void SetEntityStatsPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, const NSMesh::Ownership& p_ownership, size_t p_data );

  /*
   * Function that sets the total number of entities of a given topology stored in the entity stats class
   @param topology - Topology of the entity for which to store the data
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @param number - The data to store
  */
  void SetGlobalNumEntityPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, size_t p_data );

  /*
   * This function sets the starting and ending indexes in the entity pointer vector of the entities of a given ownership and topology
   * For instance the index of the first owned vertex and the index of the last owned vertex in the entity pointer vector of a given topology
   @param topology - Topology of the entities for which we are setting the data
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @param ownership - Ownership type for which we are setting the start and end index
   @param size - Number of entities of a given ownership used to determine the ending index in the entities vector
   @param start_index - Starting index of the entites of a given ownership type
  */
  void SetEntityOwnershipIndexesPerTopology( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName,
					     const NSMesh::Ownership& p_ownership, size_t p_size, size_t p_startIndex );

  /*
   * Add an entity to the entity vector
   @param entity - Pointer to the entity to be added
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @param topology - Topology of the entity being added
   @param gid - Global id of the new entity
   @param index - Index in the entity vector where to add the entity
  */
  void AddEntityAtIndex( NSMesh::NSEntity::Entity* p_entity, const std::string& p_interiorPartName, const NSMesh::Topology& p_topology, unsigned p_gid, size_t p_index );

  /*
   * Function that rebuilds the gid to lid map for entities of a specific topology
   @param topology - Topology of the entities whose gid to lid map is to be rebuilt
  */
  void RebuildGIDToLIDMapPerTopology( const NSMesh::Topology& p_topology );

  /*
   * Function to get the vector of pairs of start and end index in the entities vector of entities of a given ownership
   * The reason why it is a vector of pairs of index and not just a pair of index is because in the case that any new entity that is added of any ownership type will be added to the end of the entity vector
   * In this case, entites of a given owernship type will not be contiguous in the entity vector, hence, the need of a vector of pair of indexes that specifies the regions in the entity vector where entities of a given ownership type are stored
   @param topology - Topology of the entities for which we are getting the information
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @param ownership - Ownership type of the entity index information we want to get
   @output Vector where each entry is a pair of indexes corresponding to the first and last entry in the entity vector that holds entities with a given ownership type
  */
  const std::vector<std::pair<long unsigned, long unsigned> >& GetEntityIndexVec( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName,
										  const NSMesh::Ownership& p_ownership ) const;

  /*
   * Function that returns the entity vector
   @param topology - Topology of the entities for which we are getting the information   
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @output Vector of entities of a given topology
  */
  const std::vector<NSMesh::NSEntity::Entity*>& GetEntityVec( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName ) const;

  /*
   * Getter function of the number of an entity for a specified set of ownership type per processor
   @param ownership - Vector of ownership types for which to return the total number of a given entity
   @param interiorPartName - The name of the interior mesh part attached to the boundary mesh part
   @param topology - Topology of the entity we are querying
   @output Total number of a given entity of the specified onwership type per processor
  */
  unsigned GetNumberOfEntity( const NSMesh::Topology& p_topology, const std::string& p_interiorPartName, const std::vector<NSMesh::Ownership>& p_ownership ) const;

  /*
   * Getter function of the mesh element topology type
   @output Mesh element topology type
  */
  NSMesh::ElementTopology GetMeshElementTopology() const;

  /*
   * Getter function of the number of mesh vertices per element
   * This function assumes that all the mesh elements are of the same topology
   @output Number of mesh vertices per element of the mesh
  */
  int GetNumberOfVerticesPerElement() const;   

  /*
   * Setter function of the list of interior mesh parts name that are attached to this boundary mesh part
   @param attachedInteriorMeshPartsName - A vector holding the names of the interior mesh parts attached to this boundary mesh part
  */
  void SetAttachedInteriorMeshPartsName( const std::vector<std::string>& p_attachedInteriorMeshPartsName );
  
  /*
   * Getter function of the interior mesh part names attached to this boundary mesh part
   @output A vector holding the names of the interior mesh parts attached to the boundary mesh part
  */
  const std::vector<std::string>& GetAttachedInteriorMeshPartsNames() const;
  
protected:
  std::string _name = {}; /**<The name of the mesh part*/
  int _id = {}; /**<The id that defines the part, typically specified from the mesh generator or through some algorithm within the code*/
  NSMesh::Topology _primaryTopology = {}; /**<The primary topology of this part. The reason why this field is needed is mostly due to the way mesh files are written*/
  std::vector<std::map<std::string, std::unique_ptr<NSMesh::NSEntity::EntityContainer<NSMesh::NSEntity::Entity*> > > > _entityInfoPerTopology = {}; /**<A vector of a map from the name of the interior part to which the boundary is connected to an object holding information about entities of a specific topology*/
  NSMesh::ElementTopology _elementTopology = {}; /**<Topology of the elements of the subset of the mesh defined by the part*/
  int _numberOfVerticesPerElement = {}; /**<Number of mesh vertices per element. For now this is assumed to be constant for all the elements of the mesh part*/
  std::vector<std::string> _attachedInteriorMeshPartNames = {}; /**<A vector holding the names of the interior mesh parts attached to this boundary mesh part*/
};

} // namespace mesh

#endif
