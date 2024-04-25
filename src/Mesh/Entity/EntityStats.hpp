#ifndef ENTITYSTATS_HPP_
#define ENTITYSTATS_HPP_

/*
 * @brief This is a simple class that holds data about the total number of a given entity in a mesh
 */

#include<vector>
#include<map>

#include "Mesh/MeshEnum.hpp"

namespace NSMesh { namespace NSEntity {

class EntityStats {
public:
  /*
   * Default constructor
   */
  EntityStats();

  /*
   * Copy constructor not implemented
   @param eS - EntityStats to copy from
  */
  EntityStats( const EntityStats& p_eS ) = delete;

  /*
   * Equal operator not implemented
   @param eS - EntityStats to set equal to
  */
  EntityStats& operator=( const EntityStats& p_eS ) = delete;

  /*
   * Destructor
   */
  ~EntityStats();

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

private:
  std::map<NSMesh::Ownership, unsigned> _ownershipToEntityNumberMap = {}; /**<Map from the ownership type to the number of entities of a specific topology per processor*/
  size_t _globalNumberEntity = {}; /**<Total number of entities of a specific topology on the whole mesh*/
};

} // namespace Entity
} // namespace Mesh

#endif
