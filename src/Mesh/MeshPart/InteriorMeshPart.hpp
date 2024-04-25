#ifndef INTERIORMESHPART_HPP_
#define INTERIORMESHPART_HPP_

#include<string>
#include<memory>

#include "Types/Types.hpp"
#include "Mesh/Entity/EntityContainer.hpp"

namespace NSParser { namespace NSInput { class Datafile; } }

namespace NSMesh {
  
class InteriorMeshPart {
public:
  /*
   * Constructor
   */
  explicit InteriorMeshPart();

  /*
   * Copy constructor not implemented
   @param mP - Mesh part to copy from
  */
  InteriorMeshPart( const InteriorMeshPart& p_mP ) = delete;

  /*
   * Equal operator not implemented
   @param mP - Mesh part to set equal to
  */
  InteriorMeshPart& operator=( const InteriorMeshPart& p_mP ) = delete;

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - pointer to the input file
   @param section_path - a string representing the path to the section in the input file where the mesh parameters are specified
   @param section_name - name of the section in the input file that holds all the information about the mesh parameters
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

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
  NSMesh::Topology GetPrimaryTopology() const;

  /*
   * Function to allocate memory to store the entities of a topology defined on this mesh part
   @param topology - Topology of the entities for which to allocate memory
   @param size - Number of entities to allocated memory for
  */
  void AllocateMemoryPerTopology( const NSMesh::Topology& p_topology, size_t p_size );

  /*
   * Function that sets statistical data, of a given ownership, about the entities of a specific topology defined on this mesh part
   @param topology - Topology of the entities for which we are setting the statistical data
   @param ownership - Ownership type of the entities for which we are storing the information
   @param data - Number of entities of a given ownership and of a given topology
  */
  void SetEntityStatsPerTopology( const NSMesh::Topology& p_topology, const NSMesh::Ownership& p_ownership, size_t p_data );

  /*
   * Function that sets the total number of entities of a given topology stored in the entity stats class
   @param topology - Topology of the entity for which to store the data
   @param number - The data to store
  */
  void SetGlobalNumEntityPerTopology( const NSMesh::Topology& p_topology, size_t p_data );

  /*
   * This function sets the starting and ending indexes in the entity pointer vector of the entities of a given ownership and topology
   * For instance the index of the first owned vertex and the index of the last owned vertex in the entity pointer vector of a given topology
   @param topology - Topology of the entities for which we are setting the data
   @param ownership - Ownership type for which we are setting the start and end index
   @param size - Number of entities of a given ownership used to determine the ending index in the entities vector
   @param start_index - Starting index of the entites of a given ownership type
  */
  void SetEntityOwnershipIndexesPerTopology( const NSMesh::Topology& p_topology, const NSMesh::Ownership& p_ownership, size_t p_size, size_t p_startIndex );

  /*
   * Add an entity to the entity vector
   @param entity - Pointer to the entity to be added
   @param topology - Topology of the entity being added
   @param gid - Global id of the new entity
   @param index - Index in the entity vector where to add the entity
  */
  void AddEntityAtIndex( NSMesh::NSEntity::Entity* p_entity, const NSMesh::Topology& p_topology, unsigned p_gid, size_t p_index );

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
   @param ownership - Ownership type of the entity index information we want to get
   @output Vector where each entry is a pair of indexes corresponding to the first and last entry in the entity vector that holds entities with a given ownership type
  */
  const std::vector<std::pair<long unsigned, long unsigned> >& GetEntityIndexVec( const NSMesh::Topology& p_topology, const NSMesh::Ownership& p_ownership ) const;

  /*
   * Function that returns the entity vector
   @output Vector of entities of a given topology
  */
  const std::vector<NSMesh::NSEntity::Entity*>& GetEntityVec( const NSMesh::Topology& p_topology ) const;

  /*
   * Getter function of the number of an entity for a specified set of ownership type per processor
   @param ownership - Vector of ownership types for which to return the total number of a given entity
   @param topology - Topology of the entity we are querying
   @output Total number of a given entity of the specified onwership type per processor
  */
  unsigned GetNumberOfEntity( const NSMesh::Topology& p_topology, const std::vector<NSMesh::Ownership>& p_ownership ) const;

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
   * Function that sets the information about the topology of the elements that make up the part
   @param numVerticesPerElement - The number of vertices per element
   @param elementTopology - An enum that specifies the shape of the element
  */
  void SetElementTopologyInfo( int p_numVerticesPerElement, const NSMesh::ElementTopology& p_elementTopology );

  /*
   * Function that checks whether an entity belongs to the container or not
   @param gid - The global id of the entity being checked
   @param topology - The topology of the entity being checked
   @output Boolean specifying whether the entity is stored or not in the container
  */
  bool IsPresent( unsigned p_gid, const NSMesh::Topology& p_topology ) const;
  
protected:
  std::string _name = {}; /**<The name of the mesh part*/
  int _id = {}; /**<The id that defines the part, typically specified from the mesh generator or through some algorithm within the code*/
  std::vector<std::unique_ptr<NSMesh::NSEntity::EntityContainer<NSMesh::NSEntity::Entity*> > > _entityInfoPerTopology = {}; /**<One dimensional array of entity conatiners holding information about entities of a specific topology defined on this mesh part*/
  NSMesh::ElementTopology _elementTopology = {}; /**<Topology of the elements of the subset of the mesh defined by the part*/
  int _numberOfVerticesPerElement = {}; /**<Number of mesh vertices per element. For now this is assumed to be constant for all the elements of the mesh part*/
};

} // namespace mesh

#endif
