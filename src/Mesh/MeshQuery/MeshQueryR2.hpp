#ifndef MESHQUERYR2_HPP_
#define MESHQUERYR2_HPP_

/*
 * @brief Particular implementation of the mesh query class
 */

#include "Mesh/MeshQuery/MeshQuery.hpp"
#include "Factory/Factory.hpp"

namespace NSMesh {

namespace NSEntity { class Entity; }

class MeshQueryR2 : public MeshQuery {
public:
  /*
   * Default constructor
   */
  explicit MeshQueryR2();

  /*
   * Copy constructor not implemented
   @param meshQueryR2 - Mesh query to copy from
  */
  MeshQueryR2( const MeshQueryR2& p_meshQueryR2 ) = delete;

  /*
   * Equal operator not implemented
   @param meshQueryR2 - Equal operator to set equal to
  */
  MeshQueryR2& operator=( const MeshQueryR2& p_meshQueryR2 ) = delete;

  /*
   * Destructor
   */
  ~MeshQueryR2();

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - pointer to the input file
   @param section_path - a string representing the path to the section in the input file where the mesh query parameters are specified
   @param section_name - name of the section in the input file that holds all the information about the mesh query parameters
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) override;

  /*
   * Initialization function
   @param mpi_com - The communicator that will be invoked whenever any parallel communication need to be done
   @param Ndim - Dimension of the problem
  */
  void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_mpiCom, int p_nDdim ) override;

  /*
   * This is a connectivity related routine. It returns a vector of pointers to entities of a given topology that are connected to entity
   @param from_topology - Topology of the connected entities we are asking for
   @param to_topology - Topology of the connected entities we are asking for
   @param p_connectivityInfo - 2D array containing information about the connected entities of a given entity
   @param p_connectedEntitiesVec - 1D array to hold the pointers to the connected entities of a given topolgoy to some entity
  */
  void Begin( const NSMesh::Topology& p_fromTopology, const NSMesh::Topology& p_toTopology, const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const override;

private:
  /*
   * Initialization of the map that stores the connectivity info from an entity to other entities of a different topology
  */
  void _initializeConnectivityInfoProtocol();

  /*
   * Initialization of the array that holds the function pointers to the functions reponsible for the connectivity query operations
   */
  void _initializeConnectivityFunctionPointersArray();

  /*
   * This routine returns the element itslef
   @param p_connectivityInfo - Explicitly stored connectivity relations for the element we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _elementToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the faces attached to a given element
   @param p_connectivityInfo - Explicitly stored connectivity relations for the element we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _elementToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the edges attached to a given element
   @param p_connectivityInfo - Explicitly stored connectivity relations for the element we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _elementToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the vertices attached to a given element
   @param p_connectivityInfo - Explicitly stored connectivity relations for the element we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _elementToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the face itself
   @param p_connectivityInfo - Explicitly stored connectivity relations for the face we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _faceToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the elements attached to a given face
   @param p_connectivityInfo - Explicitly stored connectivity relations for the face we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _faceToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the edges attached to a given face
   @param p_connectivityInfo - Explicitly stored connectivity relations for the face we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _faceToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the vertices attached to a given face
   @param p_connectivityInfo - Explicitly stored connectivity relations for the face we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */
  void _faceToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the edge itself
   @param p_connectivityInfo - Explicitly stored connectivity relations for the edge we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _edgeToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the elements attached to a given edge
   @param p_connectivityInfo - Explicitly stored connectivity relations for the edge we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _edgeToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the faces attached to a given edge
   @param p_connectivityInfo - Explicitly stored connectivity relations for the edge we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _edgeToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the vertices attached to a given edge
   @param p_connectivityInfo - Explicitly stored connectivity relations for the edge we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _edgeToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the vertex itself
   @param p_connectivityInfo - Explicitly stored connectivity relations for the vertex we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _vertexToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the elements attached to a given vertex
   @param p_connectivityInfo - Explicitly stored connectivity relations for the vertex we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _vertexToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the faces attached to a given vertex
   @param p_connectivityInfo - Explicitly stored connectivity relations for the vertex we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _vertexToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

  /*
   * This routine returns the edges attached to a given vertex
   @param p_connectivityInfo - Explicitly stored connectivity relations for the vertex we wish to return some other connected entities to it
   @param p_connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
  */  
  void _vertexToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const;

protected:
  typedef void (MeshQueryR2::*ConnectivityFunction)( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >&, std::vector<const NSMesh::NSEntity::Entity*>& ) const; /**<Typedef to define a function pointer that takes two arguments*/
  std::vector<std::vector<ConnectivityFunction> > _getConnectivity = {}; /**<2D array holding pointers to functions responsible for running the connectivity query operations */
};

/*
* Create function used by the factor to create this class
@output A mesh query whose static type is MeshQuery and whose dynamic type is MeshQueryR2
*/
MeshQuery* createMeshQueryR2();

namespace {
  static bool register_MeshQueryR2 = MeshQuery::factory_mesh_query_Type::Register( "MESHQUERYR2", &createMeshQueryR2 );
}

} // namespace NSMesh

#endif
