#include "Mesh/MeshQuery/MeshQueryR2.hpp"
#include "Utility/Constants.hpp"
#include "Mesh/Entity/Entity.hpp"

#include<algorithm>

namespace NSMesh {

MeshQuery* createMeshQueryR2() {
  return new MeshQueryR2;
}


MeshQueryR2::MeshQueryR2() {
}

  
MeshQueryR2::~MeshQueryR2() {
}


void
MeshQueryR2::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  MeshQuery::Setup( p_datafile, p_sectionPath, p_sectionName );
}


void
MeshQueryR2::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_mpiCom, int p_nDim ) {
  MeshQuery::Initialize( p_mpiCom, p_nDim );

  // Define the connectivity protocole, i.e, the connectivity information stored between the different entities
  _initializeConnectivityInfoProtocol();

  // Initialize the 2D array of function pointers to the connectivity operations
  _initializeConnectivityFunctionPointersArray();
}


void
MeshQueryR2::_initializeConnectivityInfoProtocol() {
  _connectivityInfoProtocol.insert( std::make_pair( NSMesh::Topology::tp_ELEMENT, std::vector<NSMesh::Topology>{NSMesh::Topology::tp_VERTEX} ) );
  if( _ndim == 3 ) {
    _connectivityInfoProtocol.find( NSMesh::Topology::tp_ELEMENT )->second.push_back( NSMesh::Topology::tp_FACE );
    _connectivityInfoProtocol.insert( std::make_pair( NSMesh::Topology::tp_FACE, std::vector<NSMesh::Topology>{NSMesh::Topology::tp_ELEMENT, NSMesh::Topology::tp_EDGE} ) );
  }
  if( _ndim > 1 ) {
    _connectivityInfoProtocol.insert( std::make_pair( NSMesh::Topology::tp_EDGE, std::vector<NSMesh::Topology>{NSMesh::Topology::tp_VERTEX} ) );
  }
  _connectivityInfoProtocol.insert( std::make_pair( NSMesh::Topology::tp_VERTEX, std::vector<NSMesh::Topology>{NSMesh::Topology::tp_ELEMENT, NSMesh::Topology::tp_EDGE} ) );
}


void
MeshQueryR2::_initializeConnectivityFunctionPointersArray() {
  // The idea to build this 2D array is to avoid conditional statements since the function that will be used to retrieve the connectivity information
  // will depend on the entity topology we are interested in as well as the topology of the entities connected to it that we want to get
  // The 2D array will have the following:
  // dim1: index based on the topology of the entity of interest
  // dim2: index based on the topology of the connected entity of interest
  // output: the array will store function pointers to the corresponding functions that will return the connected entities to a given entity

  // Get the total number of mesh entities currently supported
  // The +1 is because the numbering of the enumeration class starts from 0
  int numMeshEntities = NSMesh::Topology::tp_LAST + 1;
  _getConnectivity.resize( numMeshEntities );

  // Add the element to other entities connectivity function pointers
  std::vector<ConnectivityFunction> element_connectivity_function_vector( numMeshEntities, nullptr );
  element_connectivity_function_vector[ NSMesh::Topology::tp_ELEMENT ] = &MeshQueryR2::_elementToElement;
  element_connectivity_function_vector[ NSMesh::Topology::tp_FACE ] = &MeshQueryR2::_elementToFace;
  element_connectivity_function_vector[ NSMesh::Topology::tp_EDGE ] = &MeshQueryR2::_elementToEdge;
  element_connectivity_function_vector[ NSMesh::Topology::tp_VERTEX ] = &MeshQueryR2::_elementToVertex;

  // Set the _getConnectivity entry corresponding to the element connectivity query functions
  _getConnectivity[ NSMesh::Topology::tp_ELEMENT ] = element_connectivity_function_vector;

  // Add the face to other entities connectivity function pointers
  std::vector<ConnectivityFunction> face_connectivity_function_vector( numMeshEntities, nullptr );
  face_connectivity_function_vector[ NSMesh::Topology::tp_FACE ] = &MeshQueryR2::_faceToFace;
  face_connectivity_function_vector[ NSMesh::Topology::tp_ELEMENT ] = &MeshQueryR2::_faceToElement;
  face_connectivity_function_vector[ NSMesh::Topology::tp_EDGE ] = &MeshQueryR2::_faceToEdge;
  face_connectivity_function_vector[ NSMesh::Topology::tp_VERTEX ] = &MeshQueryR2::_faceToVertex;

  // Set the _getConnectivity entry corresponding to the face connectivity query functions
  _getConnectivity[ NSMesh::Topology::tp_FACE ] = face_connectivity_function_vector;

  // Add the edge to other entities connectivity function pointers
  std::vector<ConnectivityFunction> edge_connectivity_function_vector( numMeshEntities, nullptr );
  edge_connectivity_function_vector[ NSMesh::Topology::tp_EDGE ] = &MeshQueryR2::_edgeToEdge;
  edge_connectivity_function_vector[ NSMesh::Topology::tp_ELEMENT ] = &MeshQueryR2::_edgeToElement;
  edge_connectivity_function_vector[ NSMesh::Topology::tp_FACE ] = &MeshQueryR2::_edgeToFace;
  edge_connectivity_function_vector[ NSMesh::Topology::tp_VERTEX ] = &MeshQueryR2::_edgeToVertex;

  // Set the _getConnectivity entry corresponding to the edge connectivity query functions
  _getConnectivity[ NSMesh::Topology::tp_EDGE ] = edge_connectivity_function_vector;

  // Add the vertex to other entities connectivity function pointers
  std::vector<ConnectivityFunction> vertex_connectivity_function_vector( numMeshEntities, nullptr );
  vertex_connectivity_function_vector[ NSMesh::Topology::tp_VERTEX ] = &MeshQueryR2::_vertexToVertex;
  vertex_connectivity_function_vector[ NSMesh::Topology::tp_ELEMENT ] = &MeshQueryR2::_vertexToElement;
  vertex_connectivity_function_vector[ NSMesh::Topology::tp_FACE ] = &MeshQueryR2::_vertexToFace;
  vertex_connectivity_function_vector[ NSMesh::Topology::tp_EDGE ] = &MeshQueryR2::_vertexToEdge;

  // Set the _getConnectivity entry corresponding to the vertex connectivity query functions
  _getConnectivity[ NSMesh::Topology::tp_VERTEX ] = vertex_connectivity_function_vector;
}


void
MeshQueryR2::Begin( const NSMesh::Topology& p_fromTopology, const NSMesh::Topology& p_toTopology, const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  (this->*_getConnectivity[ p_fromTopology ][ p_toTopology ])( p_connectivityInfo, p_connectedEntitiesVec );
}


void
MeshQueryR2::_elementToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_ELEMENT ];
}


void
MeshQueryR2::_elementToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_FACE ];
}


void
MeshQueryR2::_elementToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  // Get the vertices attached to the element
  const std::vector<const NSMesh::NSEntity::Entity*>& verticesVector = p_connectivityInfo[ NSMesh::Topology::tp_VERTEX ];

  // Loop over the vertices
  std::map<const NSMesh::NSEntity::Entity*, int> edge_cnt_map;
  for( size_t n=0; n<verticesVector.size(); ++n ) {
    // For each vertex
    const NSMesh::NSEntity::Entity& vertex = *verticesVector[n];

    // Get the edges attached to the vertex
    if( vertex.GetNumConnectivity( NSMesh::Topology::tp_EDGE ) > 0 ) {
      const auto& edgesVector = vertex._connectivity[ NSMesh::Topology::tp_EDGE ];
      
      // Loop over the egdes
      for( size_t e=0; e<edgesVector.size(); ++e ) {
	if( edge_cnt_map.find( edgesVector[e] ) == edge_cnt_map.end() ) {
	  // If it is the first time this edge is encountered add to the map
	  edge_cnt_map.insert( std::make_pair( edgesVector[e], 1 ) );
	}
	else {
	  // If it is not the first time the edge is encountered increment its counter
	  edge_cnt_map.find( edgesVector[e] )->second++;
	}
      }
    }
  }

  // Return the edges in the map whose counter is greater than 1
  std::map<const NSMesh::NSEntity::Entity*, int>::const_iterator it = edge_cnt_map.begin();
  while( it != edge_cnt_map.end() ) {
    if( it->second > 1 ) {
      p_connectedEntitiesVec.push_back( it->first );
    }

    it++;
  }
}


void
MeshQueryR2::_elementToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_VERTEX ];
}


void
MeshQueryR2::_faceToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_FACE ];
}


void
MeshQueryR2::_faceToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_ELEMENT ];
}


void
MeshQueryR2::_faceToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_EDGE ];
}


void
MeshQueryR2::_faceToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  // Get the edges forming the face
  const std::vector<const NSMesh::NSEntity::Entity*>& edges_vector = p_connectivityInfo[ NSMesh::Topology::tp_EDGE ];

  // Loop over the edges
  for( size_t e=0; e<edges_vector.size(); ++e ) {
    // For each edge
    const NSMesh::NSEntity::Entity& edge = *edges_vector[e];

    // Get the two vertices forming the edge
    const auto& vertices_vector = edge._connectivity[ NSMesh::Topology::tp_VERTEX ];

    for( size_t n=0; n<vertices_vector.size(); ++n ) {
      // Push back the vertices to the connected entities vector
      p_connectedEntitiesVec.push_back( vertices_vector[n] );
    }
  }

  // Remove the duplicated vertices in the connected entities vector
  std::sort( p_connectedEntitiesVec.begin(), p_connectedEntitiesVec.end() );
  p_connectedEntitiesVec.erase( std::unique( p_connectedEntitiesVec.begin(), p_connectedEntitiesVec.end() ), p_connectedEntitiesVec.end() );
}
 

void
MeshQueryR2::_edgeToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_EDGE ];
}

 
void
MeshQueryR2::_edgeToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  // Get the vertices forming the edge
  const std::vector<const NSMesh::NSEntity::Entity*>& verticesVector = p_connectivityInfo[ NSMesh::Topology::tp_VERTEX ];
  
  // We assume an edge has two vertices
  // For each vertex get the connected elements
  const NSMesh::NSEntity::Entity& vertex1 = *verticesVector[0];
  auto vertex1ElementsVector = vertex1._connectivity[ NSMesh::Topology::tp_ELEMENT ];

  const NSMesh::NSEntity::Entity& vertex2 = *verticesVector[1];
  auto vertex2ElementsVector = vertex2._connectivity[ NSMesh::Topology::tp_ELEMENT ];
  
  // Take the intersection of the two elements vectors
  // The result of the intersection will be the elements connected to the edge
  std::sort( vertex1ElementsVector.begin(), vertex1ElementsVector.end() );
  std::sort( vertex2ElementsVector.begin(), vertex2ElementsVector.end() );

  std::set_intersection( vertex1ElementsVector.begin(), vertex1ElementsVector.end(), vertex2ElementsVector.begin(), vertex2ElementsVector.end(), std::back_inserter( p_connectedEntitiesVec ) );
}


void
MeshQueryR2::_edgeToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  // Get the elements connected to the edge
  std::vector<const NSMesh::NSEntity::Entity*> elementsVector;
  _edgeToElement( p_connectivityInfo, elementsVector );

  // Loop over each element
  std::map<const NSMesh::NSEntity::Entity*, int> face_cnt_map;
  for( size_t e=0; e<elementsVector.size(); ++e ) {
    if( elementsVector[e]->GetNumConnectivity( NSMesh::Topology::tp_FACE ) > 0 ) {
      // Get the faces of this element
      const std::vector<const NSMesh::NSEntity::Entity*>& facesVector = elementsVector[e]->_connectivity[ NSMesh::Topology::tp_FACE ];
      
      // Loop over the faces and insert them in a map from the face to a counter representing the number of occurences
      for( size_t f=0; f<facesVector.size(); ++f ) {
	if( face_cnt_map.find( facesVector[f] ) == face_cnt_map.end() ) {
	  // If it is the first time this face is encountered add to the map
	  face_cnt_map.insert( std::make_pair( facesVector[f], 1 ) );
	}
	else {
	  // If it is not the first time the face is encountered increment its counter
	  face_cnt_map.find( facesVector[f] )->second++;
	}
      }
    }
  }

  // Return the faces in the map whose counter is greater than 1
  std::map<const NSMesh::NSEntity::Entity*, int>::const_iterator it = face_cnt_map.begin();
  while( it != face_cnt_map.end() ) {
    if( it->second > 1 ) {
      p_connectedEntitiesVec.push_back( it->first );
    }

    it++;
  }
}


void
MeshQueryR2::_edgeToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_VERTEX ];
}


void
MeshQueryR2::_vertexToVertex( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_VERTEX ];
}


void
MeshQueryR2::_vertexToElement( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_ELEMENT ];
}


void
MeshQueryR2::_vertexToFace( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  // Get the vector of elements connected to the vertex
  const std::vector<const NSMesh::NSEntity::Entity*>& elementsVector = p_connectivityInfo[ NSMesh::Topology::tp_ELEMENT ];

  // Loop over each element
  std::map<const NSMesh::NSEntity::Entity*, int> face_cnt_map;
  for( size_t e=0; e<elementsVector.size(); ++e ) {
    if( elementsVector[e]->GetNumConnectivity( NSMesh::Topology::tp_FACE ) > 0 ) {
      // Get the faces of this element
      const std::vector<const NSMesh::NSEntity::Entity*>& facesVector = elementsVector[e]->_connectivity[ NSMesh::Topology::tp_FACE ];
      
      // Loop over the faces and insert them in a map from the face to a counter representing the number of occurences
      for( size_t f=0; f<facesVector.size(); ++f ) {
	if( face_cnt_map.find( facesVector[f] ) == face_cnt_map.end() ) {
	  // If it is the first time this face is encountered add to the map
	  face_cnt_map.insert( std::make_pair( facesVector[f], 1 ) );
	}
	else {
	  // If it is not the first time the face is encountered increment its counter
	  face_cnt_map.find( facesVector[f] )->second++;
	}
      }
    }
  }

  // Return the faces in the map whose counter is greater than 1
  std::map<const NSMesh::NSEntity::Entity*, int>::const_iterator it = face_cnt_map.begin();
  while( it != face_cnt_map.end() ) {
    if( it->second > 1 ) {
      p_connectedEntitiesVec.push_back( it->first );
    }

    it++;
  }
}


void
MeshQueryR2::_vertexToEdge( const std::vector<std::vector<const NSMesh::NSEntity::Entity*> >& p_connectivityInfo, std::vector<const NSMesh::NSEntity::Entity*>& p_connectedEntitiesVec ) const {
  p_connectedEntitiesVec = p_connectivityInfo[ NSMesh::Topology::tp_EDGE ];
}

}
