#ifndef ENTITY_HPP_
#define ENTITY_HPP_

/*
 * @brief The entity class is an abstract representation of the different components that constitute a mesh
 * Different entity types currently supported are vertex, elements, edges, and faces
 * An entity is uniquely determined by its global id, topology, and its map of connectivity information
 * The global id is typically read from a mesh file or assigned through some algorithm in the case of mesh refinement for example
 * The topology sets the type of the entity (vertex, element, face or edge)
 * The connectivity map stores information about the entities of different topology that are connected to this entity
 * The amount of information stored in this map is determined in the code from the mesh query class
 */

#include<memory>
#include<vector>
#include<map>
#include<algorithm>

#include "Mesh/MeshEnum.hpp"

namespace NSMesh { class MeshQueryR0; }
namespace NSMesh { class MeshQueryR2; }

namespace NSMesh { namespace NSEntity {
    
class Entity {
public:
  /*
   * Constructor
   @param topology - The topology of the entity
   @param gid - The global id of the entity which is unique across all other entities of the mesh and all other processors
   @param lid - The local id of the entity with respect to the processor who has access to it (i.e. owns or shares)
  */
  explicit Entity( const NSMesh::Topology& p_topology, unsigned p_gid, unsigned p_lid );

  /*
   * Copy constructor
   @param e - Entity to be copied
  */
  explicit Entity( const Entity& p_e );

  /*
   * Default constructor not implemented
   */
  Entity() = delete;
  
  /*
   * Equal operator not implemented
   */
  Entity& operator=( const Entity& p_e ) = delete;

  /*
   * Function to set the global id of the entity
   @param gid - Global id
  */
  void SetGid( unsigned p_gid );

  /*
   * Function to set the local id of the entity
   @param lid - Local id
  */
  void SetLid( unsigned p_lid );

  /*
   * Function to set the topology of the entity
   @param topology - Topology to be set (vertex, element, edge, face)
  */
  void SetTopology( const NSMesh::Topology& p_topology );
  
  /*
   * Function to get the global id of the entity
   @output Global id of the entity
  */
  unsigned GetGid() const;

  /*
   * Getter function of the entity lid
   @output Local id of the entity
  */
  const unsigned& GetLid() const;
  
  /*
   * Function to get the topology of the entity
   @output Topology of the entity
  */
  NSMesh::Topology GetTopology() const;

  /*
   * Function to add a certain connectivity information about the entity
   @param topology - Topology of the connected entity that we would like to store the information about
   @param p_entity - Pointer to the connected entity that we are are storing the connectivity information about
  */
  void AddConnectivity( const NSMesh::Topology& p_topology, const NSMesh::NSEntity::Entity* p_entity );

  /*
   * Function that returns the number of connected entities of a specified topology
   @param topology - Topology of the connected entities we want to recover
   @output Number of connected entities of the specified topology
  */
  int GetNumConnectivity( const NSMesh::Topology& p_topology ) const;

  /*
   * Function that returns the connectivity information
   @output Connectivity 2D array containing the information about the connected entities to this entity
  */
  const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& GetConnectivityInfo() const;

  /*
   * Declaring MeshQuery as a friend class of Entity to be able to access the connectivity information stored in Entity
   */
  friend class NSMesh::MeshQueryR0;
  friend class NSMesh::MeshQueryR2;
  
private:
  NSMesh::Topology _topology = {}; /**<Topology of the entity (vertex, element, face, edge)*/
  unsigned _gid = {}; /**<Global id of the entity*/
  unsigned _lid = {}; /**<Local id of the entity*/
  std::vector<std::vector<const NSMesh::NSEntity::Entity*>> _connectivity = {}; /**<Two dimensional array holding information about the connected entities to this entity*/
};

} // namespace Entity
} // namespace Mesh

#endif
