#ifndef MESHQUERY_HPP_
#define MESHQUERY_HPP_

/*
 * @brief This class is created to allow storage of different information about the mesh
 * When storing a mesh, there is a tradeoff between the amount of memory consumed to store the mesh and the speed to query the mesh
 * Query functions on the mesh can be as simple as recovering connected elements to a vertex for example or adding new entities to the mesh
 * The memory and computation cost of those query functions on the mesh depends on the explicity connectivity information stored about the mesh
 * This class is an abstract class that represents the different models of mesh storage and corresponding query information
 * The connectivity query class and the naming used in this code is based on the following paper:
 * "Mesh data structure selection for mesh generation and FEA application", Rao Garimella
 */

#include<memory>
#include<map>
#include<vector>
#include<string>

#include "Mesh/MeshEnum.hpp"

namespace NSFactory { template<typename Product, typename Identifier, typename... Argument> class Factory; }
namespace NSParser { namespace NSInput { class Datafile; } }
namespace NSCommunicator { namespace NSMPIComWrapper { class MPIComWrapper; } }

namespace NSMesh {

namespace NSEntity { class Entity; }

class MeshQuery {
public:
  typedef NSFactory::Factory<MeshQuery*, std::string, void> factory_mesh_query_Type; /* *< Mesh query factor */

  /*
   * Default constructor
   */
  explicit MeshQuery();

  /*
   * Copy constructor not implemented
   @param meshQuery - Mesh query to copy from
  */
  MeshQuery( const MeshQuery& p_meshQuery ) = delete;

  /*
   * Equal operator
   @param meshQuery - Mesh query to set equal to
  */
  MeshQuery& operator=( const MeshQuery& p_meshQuery ) = delete;

  /*
   * Destructor
   */
  virtual ~MeshQuery();

  /*
   * Setup function that mainly reads the input parameters corresponding to the mesh query class from the input file
   @param datafile - pointer to the input file
   @param section_path - a string representing the path to the section in the input file where the mesh query parameters are specified
   @param section_name - name of the section in the input file that holds all the information about the mesh query parameters
  */
  virtual void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Initialization function
   @param mpi_com - The communicator that will be invoked whenever any parallel communication need to be done
   @Ndim - Dimension of the problem we are solving
  */
  virtual void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_mpiCom, int p_nDim );

  /*
   * Function that returns the information about the connectivity information to be explicitly stored by the mesh
   @ouput Map specifying for each entity with topology what are the connected entities with different topologies for which we will explicitly store the connectivity relations
  */
  const std::map<NSMesh::Topology, std::vector<NSMesh::Topology> >& GetConnectivityInfoProtocol() const;

  /*
   * This is a connectivity related routine. It returns a vector of pointers to entities of a given topology that are connected to entity
   @param topology - Topology of the connected entities we are asking for
   @param connectedEntitiesVec - Vector of pointers to the entities of the given topology that are connected to entity
   @param p_connectivityInfo - 2D array containing information about the connected entities of a given entity
   @output p_connectedEntitiesVec - 1D array to hold the pointers to the connected entities of a given topolgoy to some entity
  */
  virtual void Begin( const NSMesh::Topology& p_fromTopology, const NSMesh::Topology& p_toTopology, const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityMap, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const = 0;

  /*
   * Getter function of the spatial dimension of the problem
   @output Spatial dimension of the problem
  */
  int GetSpatialDimension() const;

  /*
   * Function that returns the type of the mesh query
   @output Type of the mesh query
  */
  NSMesh::MeshQueryType GetType() const;

protected:
  std::string _name = {}; /**<Mesh query name*/
  std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper> _mpiCom = {}; /**<MPI communicator responsible for doing any communication on the mesh*/
  int _ndim = {}; /**<Dimension of the problem we are solving*/
  std::map<NSMesh::Topology, std::vector<NSMesh::Topology> > _connectivityInfoProtocol = {}; /**<Map specifying for each entity with topology what are the connected entities with different topologies for which we will explicitly store the connectivity relations*/
  NSMesh::MeshQueryType _type = {}; /**<Type of the mesh query which dictates what connectivity information is stored*/
};

} // namespace NSMesh

#endif
