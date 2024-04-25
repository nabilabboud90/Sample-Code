#include "FESolutionToFieldDataHelper.hpp"
#include "Utility/Exceptions.hpp"
#include "Mesh/Mesh.hpp"
#include "Mesh/Entity/Bucket.hpp"
#include "Field/Field.hpp"
#include "LinearSystem/LinearSystem.hpp"
#include "Field/FunctionSpace.hpp"
#include "Field/FunctionSpaceHelper.hpp"

namespace NSSimulation { namespace NSProblem { namespace NSFE {

void SetSolutionVectorFromFieldsData( const NSMesh::Mesh* p_mesh, const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
				      const std::map<std::string, std::unique_ptr<NSField::Field<int> > >& p_nameToFieldDofIdsMap,
				      NSLinearSystem::LinearSystem* const p_linearSystem ) {
  // Loop over the unknown dof ids fields
  for( const auto& unknown_dof_ids_info : p_nameToFieldDofIdsMap ) {
    // Get the unknown dof ids field
    const auto& unknown_dof_ids_field = unknown_dof_ids_info.second;

    // Get the field
    const auto& unknown_field = p_nameToFieldMap.find( unknown_dof_ids_info.first )->second;
    
    // Get the parts on which the field is defined
    const auto& parts_names = unknown_dof_ids_field->GetPartsName();

    // Loop over each part
    for( const auto& part_name : parts_names ) {
      // Get the type of the entities on which nodes are defined (i.e., vertex, edge, face, element)
      auto unknown_field_element_type = unknown_dof_ids_field->GetFunctionSpace( part_name )->GetFieldElement();
      const auto& topologies_with_nodes = NSField::factory_get_entity_topology_with_nodes::Create( unknown_field_element_type )();

      // Loop over each entity type on which nodes are defined
      for( const auto& topology_with_nodes : topologies_with_nodes ) {
	// Get the vector of buckets holding the entities of the specific topology on which nodes are defined
	std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> > buckets;
	p_mesh->GetBucketsByInteriorPart( buckets, {NSMesh::Ownership::o_OWNED, NSMesh::Ownership::o_SHARED}, topology_with_nodes, part_name );

	// Create the object to hold the ids of the degrees of freedom defined on the entities in the bucket
	auto unknown_num_nodes_per_entity = NSField::factory_get_num_nodes_per_entity::Create( unknown_field_element_type )( topology_with_nodes );
	auto unknown_dofs_per_node = unknown_dof_ids_field->GetNumDofsPerNode();
	Kokkos::View<int***> unknown_dof_ids_per_bucket( "unknown_dof_ids_per_bucket", BUCKETSIZE, unknown_num_nodes_per_entity, unknown_dofs_per_node );
	Kokkos::View<double***> unknown_values_per_bucket( "unknown_values_per_bucket", BUCKETSIZE, unknown_num_nodes_per_entity, unknown_dofs_per_node );
	
	// Loop over each bucket
	for( const auto& bucket : buckets ) {
	  // Gather the ids of the dofs for the current unknown defined on the entities in the current bucket
	  unknown_dof_ids_field->GatherData( bucket, unknown_dof_ids_field->GetFunctionSpace( part_name ), NSField::AccessPer::ad_NODE, unknown_dof_ids_per_bucket );
	  
	  // Gather the unknown field values for the dofs in the current bucket
	  unknown_field->GatherData( bucket, unknown_field->GetFunctionSpace( part_name ), NSField::AccessPer::ad_NODE, unknown_values_per_bucket );

	  // Loop over each entity
	  auto num_entities = bucket.GetSize();
	  std::vector<int> dof_ids( num_entities * unknown_num_nodes_per_entity * unknown_dofs_per_node, 0 );
	  std::vector<double> dof_values( num_entities * unknown_num_nodes_per_entity * unknown_dofs_per_node, 0.0 );
	  for( size_t e=0; e<num_entities; ++ e ) {
	    for( unsigned n=0; n<unknown_num_nodes_per_entity; ++n ) {
	      for( int d=0; d<unknown_dofs_per_node; ++d ) {
		dof_ids[ e * unknown_num_nodes_per_entity * unknown_dofs_per_node + n * unknown_dofs_per_node + d ] = unknown_dof_ids_per_bucket( e, n , d );
		dof_values[ e * unknown_num_nodes_per_entity * unknown_dofs_per_node + n * unknown_dofs_per_node + d ] = unknown_values_per_bucket( e, n , d );
	      }
	    }
	  }

	  // Set the values in the solution vector of the linear system
	  p_linearSystem->SetVector( 1, dof_ids.size(), dof_ids.data(), dof_values.data() );
	}
      }
    }
  }
}


void UpdateFieldDataFromSolutionIncrement( const NSMesh::Mesh* const p_mesh, const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
					   const std::map<std::string, std::unique_ptr<NSField::Field<int> > >& p_nameToFieldDofIdsMap,
					   NSLinearSystem::LinearSystem* const p_linearSystem ) {
  // Get the solution increment field
  auto& solution_increment = *p_linearSystem->GetVector( 2 );
  
  // Loop over the unknown dof ids fields
  for( const auto& unknown_dof_ids_info : p_nameToFieldDofIdsMap ) {
    // Get the unknown dof ids field
    const auto& unknown_dof_ids_field = unknown_dof_ids_info.second;

    // Get the field
    const auto& unknown_field = p_nameToFieldMap.find( unknown_dof_ids_info.first )->second;
    
    // Get the mesh parts on which is defined in the equations
    // Note this list can be the same as that as the field is defined on or a subset of it
    const auto& parts_names = unknown_dof_ids_info.second->GetPartsName();
    
    // Build a map from the function space to the list of part names on which it is defined
    std::map<std::string, std::vector<std::string> > function_space_to_parts;
    std::map<std::string, const NSField::FunctionSpace* > name_to_fs;
    for( const auto& part_name : parts_names ) {
      // Get the pointer to the function space
      const auto& function_space = unknown_dof_ids_field->GetFunctionSpace( part_name );
      
      // Get the name of the function space
      const auto& name = function_space->GetName();

      // Update the map
      if( function_space_to_parts.find( name ) == function_space_to_parts.end() ) {
	function_space_to_parts.insert( std::make_pair( name, std::vector<std::string>( 1, part_name ) ) );
	name_to_fs.insert( std::make_pair( name, function_space ) );
      }
      else {
	function_space_to_parts.find( name )->second.push_back( part_name );
      }
    }

    // Loop over each function space
    for( const auto& info : function_space_to_parts ) {
      // Get the parts on which this function space is defined
      const auto& part_names = info.second;

      // Get the pointer to the function space
      const auto& function_space = name_to_fs.find( info.first )->second;

      // Get the field element type on the current part
      const auto& field_element_type = function_space->GetFieldElement();
      
      // Get the type of the entities on which nodes are defined (i.e., vertex, edge, face, element)
      const auto& topologies_with_nodes = NSField::factory_get_entity_topology_with_nodes::Create( field_element_type )();

      // Loop over each entity type on which nodes are defined
      for( const auto& topology_with_nodes : topologies_with_nodes ) {
	// Get some information about the field element of the current function space
	auto unknown_num_nodes_per_entity = NSField::factory_get_num_nodes_per_entity::Create( field_element_type )( topology_with_nodes );
	auto unknown_dofs_per_node = unknown_dof_ids_field->GetNumDofsPerNode();
	
	// Create the object to hold the ids of the degrees of freedom defined on the entities in the bucket
	Kokkos::View<int***> unknown_dof_ids_per_bucket( "unknown_dof_ids_per_bucket", BUCKETSIZE, unknown_num_nodes_per_entity, unknown_dofs_per_node );
	Kokkos::View<double***> unknown_values_per_bucket( "unknown_values_per_bucket", BUCKETSIZE, unknown_num_nodes_per_entity, unknown_dofs_per_node );
	
	// Loop over each part
	std::unordered_map<unsigned, std::vector<double> > entity_to_data_values;
	for( const auto& part_name : part_names ) {
	  std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> > buckets;
	  p_mesh->GetBucketsByInteriorPart( buckets, {NSMesh::Ownership::o_OWNED, NSMesh::Ownership::o_SHARED}, topology_with_nodes, part_name );
	  
	  // Loop over the buckets
	  for( size_t b=0; b<buckets.size(); ++b ) {
	    // Get the bucket
	    const auto& bucket = buckets[b];
	    
	    // Gather the ids of the dofs of the nodes defined on the entities in this bucket
	    unknown_dof_ids_field->GatherData( bucket, unknown_dof_ids_field->GetFunctionSpace( part_name ), NSField::AccessPer::ad_NODE, unknown_dof_ids_per_bucket );

	    // Gather the unknown field data at the nodes defined on the entities in this bucket
	    unknown_field->GatherData( bucket, unknown_field->GetFunctionSpace( part_name ), NSField::AccessPer::ad_NODE, unknown_values_per_bucket );

	    // Loop over each node and update the of the dofs of each node
	    for( size_t e=0; e<bucket.GetSize(); ++e ) {
	      // Get the entity
	      const auto& entity = *bucket[e];

	      // Get the entity gid
	      const auto& entity_gid = entity.GetGid();
	      
	      // Allocate memory for the data of the nodes defiend on the entity
	      entity_to_data_values.insert( std::make_pair( entity_gid, std::vector<double>( unknown_num_nodes_per_entity * unknown_dofs_per_node, 0.0 ) ) );

	      // Get the vector holding the data for the current entity
	      auto& entity_data = entity_to_data_values.find( entity_gid )->second;
	      for( unsigned n=0; n<unknown_num_nodes_per_entity; ++n ) {
		for( int d=0; d<unknown_dofs_per_node; ++d ) {
		  entity_data[ n * unknown_dofs_per_node + d ] = unknown_values_per_bucket( e, n, d ) + solution_increment[ unknown_dof_ids_per_bucket( e, n, d ) ];
		}
	      }
	    }
	  }

	  // Create a vector with the global ids of the entities on which the nodes of the function space are defined
	  std::vector<unsigned> entities_gids( entity_to_data_values.size(), 0 );
	  unsigned cnt = 0;
	  for( const auto& info : entity_to_data_values ) {
	    entities_gids[ cnt ] = info.first;
	    cnt++;
	  }

	  // Get a vector of buckets holding the entities whose global ids are contained in entities_gids
	  std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> > buckets_from_gids;
	  p_mesh->GetBucketsByGid( buckets_from_gids, entities_gids, topology_with_nodes );
	  
	  // Loop over each bucket
	  for( size_t b=0; b<buckets_from_gids.size(); ++b ) {
	    // Get the bucket
	    const auto& bucket = buckets_from_gids[b];

	    // Copy the values of the dofs of the nodes defined on the entities in the bucket
	    for( size_t e=0; e<bucket.GetSize(); ++e ) {
	      // Get the entity
	      const auto& entity = *bucket[e];

	      // Get the gid of the entity
	      const auto& entity_gid = entity.GetGid();

	      // Get the data vector of the entity
	      const auto& entity_data = entity_to_data_values.find( entity_gid )->second;

	      // Copy the data corresponding to the entity
	      for( unsigned n=0; n<unknown_num_nodes_per_entity; ++n ) {
		for( int d=0; d<unknown_dofs_per_node; ++d ) {
		  unknown_values_per_bucket( e, n, d ) = entity_data[ n * unknown_dofs_per_node + d ];
		}
	      }
	    }

	    // Scatter the data to the unknown field
	    unknown_field->ScatterData( bucket, unknown_field->GetFunctionSpace( part_name ), NSField::AccessPer::ad_NODE, unknown_values_per_bucket );
	  }
	}
      }
    }

    // Parallel communicate the data between the processors
    unknown_field->CommunicateInterfaceData( {NSMesh::Ownership::o_SHARED}, NSField::CombineType::c_REPLACE );
  }
}    
      
} // namespace NSFE
} // namespace NSProblem
} // namespace NSSimulation
