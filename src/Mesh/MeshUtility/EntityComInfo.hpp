#ifndef ENTITYCOMINFO_HPP_
#define ENTITYCOMINFO_HPP_

#include<map>
#include<vector>

#include "Types/Types.hpp"

namespace NSMesh { namespace NSFields {

class EntityComInfo {
public:
  /*
   * Default constructor
   */
  EntityComInfo();

  /*
   * Copy constructor not implemented
   */
  EntityComInfo( const EntityComInfo& p_entityComInfo );

  /*
   * Equal operator not implemented
   */
  EntityComInfo& operator=( const EntityComInfo& p_entityComInfo );
  
  /*
   * Destructor
   */
  ~EntityComInfo();

  /*
   * Getter function of ID of the owning processor of this entity
   */
  const int& GetOwningProcID() const;

  /*
   * Getter function of ID of the sharing processors of this entity
   */
  const std::vector<int>& GetSharingProcID() const;

  /*
   * Getter function of ID of the universal processors of this entity
   */
  const std::vector<int>& GetUniversalProcID() const;

  /*
   * Setter function of the ID of the owning processor
   @param p_id - Id of the owning processor
  */
  void SetOwningProcID( int p_id );

  /*
   * Setter function of the IDs of the sharing processors
   @param p_ids - Ids of the sharing processors
  */
  void SetSharingProcID( const std::vector<int>& p_ids );

  /*
   * Setter function of the IDs of the universal processors
   @param p_ids - Ids of the universal processors
  */
  void SetUniversalProcID( const std::vector<int>& p_ids );

  /*
   * Add a processor ID to the already existing vector of the sharing processors id
   @param p_id - The new id to be added
  */
  void AddSharingProcID( int p_id );

  /*
   * Add a processor ID to the already existing vector of the universal processors id
   @param p_id - The new id to be added
  */
  void AddUniversalProcID( int p_id );
  
private:
  int _ownerProc = {}; /**<ID of the owning processor of the entityt*/
  std::vector<int> _sharingProcs = {}; /**<IDs of the sharing processors of the entity*/
  std::vector<int> _universalProcs = {}; /**<IDs of the universal processors of the entity*/
};

} // namespace MeshField
} // namespace Mesh

#endif
