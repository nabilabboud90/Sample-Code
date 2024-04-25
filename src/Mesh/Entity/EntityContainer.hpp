#ifndef ENTITYCONTAINER_HPP_
#define ENTITYCONTAINER_HPP_

/*
 * @brief EntityContainer is class responsible for storing the data about all the entities of a given topology (node, element, face or edge)
 * An entity container cannot hold at the same time data about an element and a node for example
 * An entity container holds the following information
 * A topology that specifies the type of entities being references by this entity container
 * An entity vector that holds all the entities of a given topology
 * The entities in this vector are stored in a particular order as follow
 * First all the owned entities of a given topology are stored (for instance all the owned vertices)
 * Second all the shared entities of a given topology are stored (for instance all the shared vertices)
 * Third all the universsal entities of a given topology are stored (for instance all the universal vertex)
 * To keep track of the indexes of the entity vector where entities of a given ownership type ends and entities of a different ownership type starts an additional information is stored in the entity container
 * This additional information is stored in the ownership to entity indexes map which maps from an ownership type to a vector of pairs of long unsigned int where those pair of long unsigned int represent the start and end index respectively of those entities in the entity vector
 * The reason i use a vector of pair of long unsigned and not just a pair of long unsigned int is because if we later decide to add entities to the mesh (for instance mesh refinement) then we would like to add those entities to the end of the entities vector in order to avoid memory reallocation as much as possible, hence, in this case the entities of a given ownership type won't be contiguous anymore in the entity vector
 */

#include<map>
#include<unordered_map>
#include<vector>
#include<memory>

#include "Mesh/MeshEnum.hpp"

namespace NSMesh { namespace NSEntity {

class Entity;
class EntityStats;

template<typename T>
class EntityContainer {
public:
  /*
   * Default constructor not implemented
   */
  EntityContainer() = delete;

  /*
   * Constructor
   @param topology - Topology of the entities to be referenced by this entity container
  */
  explicit EntityContainer( const NSMesh::Topology& p_topology );

  /*
   * Copy constructor not implemented
   @param eC - Entity container to copy from
  */
  EntityContainer( const EntityContainer& p_eC ) = delete;

  /*
   * Equal operator not implemented
   @param eC - Entity container to set equal to
  */
  EntityContainer& operator=( const EntityContainer& p_eC ) = delete;

  /*
   * Destructor
   */
  ~EntityContainer();

  /*
   * @brief The way the implementation of this class works is to allocate first all the memory needed to store the entities in this container
   * Once the memory is allocated, loop over the entities and add them one by one
   */

  /*
   * Function to allocate memory for a given number of entities
   @param size - Number of entities to allocate memory for
  */
  void AllocateSizeForEntity( size_t p_size );

  /*
   * This function is to set the starting and ending index in the entity vec of the entites of a given ownership type
   * For instance the index of the first owned vertex and the index of the last owned vertex in the entity vector
   @param ownership - Ownership type for which we are setting the start and end index
   @param size - Number of entities of a given ownership used to determine the ending index in the entities vector
   @param start_index - Starting index of the entites of a given ownership type
  */
  void SetEntityOwnershipIndexes( const NSMesh::Ownership& p_ownership, size_t p_size, size_t p_startIndex );

  /*
   * Add an entity to the entity vector
   @param entity - A reference to the entity to be added
   @param gid - Global id of the new entity
   @param index - Index in the entity vector where to add the entity
  */
  void AddEntityAtIndex( T& p_entity, unsigned p_gid, size_t p_index );

  /*
   * Function to get the vector of pairs of start and end index in the entities vector of entities of a given ownership
   * The reason why it is a vector of pairs of index and not just a pair of index is because in the case that any new entity that is added of any ownership type will be added to the end of the entity vector
   * In this case, entites of a given owernship type will not be contiguous in the entity vector, hence, the need of a vector of pair of indexes that specifies the regions in the entity vector where entities of a given ownership type are stored
   @param ownership - Ownership type of the entity index information we want to get
   @output Vector where each entry is a pair of indexes corresponding to the first and last entry in the entity vector that holds entities with a given ownership type
  */
  const std::vector<std::pair<long unsigned, long unsigned> >& GetEntityIndexVec( const NSMesh::Ownership& p_ownership ) const;

  /*
   * Function that returns the entity vector
   @output Vector of entities of a given topology
  */
  const std::vector<T>& GetEntityVec() const;

  /*
   * Function that returns an entity given a global id
   @param gid - Globl id of the entity we are asking for
   @output The entity we asked for through the gid
  */
  const T& GetEntity( const unsigned& p_gid ) const;

  /*
   * Getter function of the number of an entity for a specified set of ownership type per processor
   @param ownership - Vector of ownership types for which to return the total number of a given entity
   @output Total number of a given entity of the specified onwership type per processor
  */
  size_t GetNumberOfEntity( const std::vector<NSMesh::Ownership>& p_ownership ) const;

  /*
   * Setter function of the number of an entity for a specified set of ownership type per processor
   @param ownership - Vector of ownership types for which the total number of a given entity is set
   @param Total number of a given entity of the specified onwership type per processor
  */
  void SetNumberOfEntity( const NSMesh::Ownership& p_ownership, size_t p_number );

  /*
   * Getter function of the total number of an entity on the whole mesh
   @output Total number of a given entity on the whole mesh
  */
  const size_t& GetGlobalNumberOfEntity() const;

  /*
   * Setter function of the total number of an entity on the whole mesh
   @param number - Total number of a given entity on the whole mesh
  */
  void SetGlobalNumberOfEntity( const size_t& p_number );

  /*
   * Function that rebuilds the gid to lid map
   */
  void RebuildGIDToLIDMap();

  /*
   * Function that checks whether an entity belongs to the container or not
   @param gid - The global id of the entity being checked
   @output Boolean specifying whether the entity is stored or not in the container
  */
  bool IsPresent( unsigned p_gid ) const;

  const std::unordered_map<unsigned, unsigned>& GetGidToLidMap() const { return _gidToLidMap; }
  
private:
  NSMesh::Topology _topology = {}; /**<Topology of the entities stored in this entity container*/
  std::map<NSMesh::Ownership, std::vector<std::pair<long unsigned, long unsigned> > > _ownershipToEntityIndexesMap = {}; /**<Map from a given ownership type to the starting and ending index in the entity vector of the entities of this ownership type. More information is given in the description at the begining of the file*/
  std::vector<T> _entityVec = {}; /**<Vector that holds all the entities (owned, shared, universal) of a given topology*/
  std::unordered_map<unsigned, unsigned> _gidToLidMap = {}; /**<Map from the global id of the entity to its local id which corresponds to its index in the entity vector. I chose to use unordered_map becasue the lookup operation is much faster compared to an std::map even though it consumes more memory (hash table vs binary search tree)*/
  std::unique_ptr<NSMesh::NSEntity::EntityStats> _entityStats; /**<Pointer to a class holding data about entities of a specific topology*/
};

} // namespace Entity
} // namespace Mesh

#endif
