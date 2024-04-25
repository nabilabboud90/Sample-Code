#include "FEPrescribedConditionsHelper.hpp"
#include "PrescribedCondition/PrescribedCondition.hpp"
#include "Utility/Exceptions.hpp"
#include "Mesh/Mesh.hpp"
#include "Mesh/Entity/Bucket.hpp"
#include "Field/Field.hpp"
#include "InputFunction/InputFunction.hpp"
#include "LinearSystem/LinearSystem.hpp"
#include "Field/FunctionSpaceHelper.hpp"
#include "Field/FunctionSpace.hpp"

namespace NSSimulation { namespace NSProblem { namespace NSFE {

void SetValuesPerBucket( const NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*>& p_bucket, unsigned p_fieldNumNodesPerElement, int p_coordinatesNumDofsPerNode,
			 int p_fieldNumDofsPerNode, std::vector<double>& p_functionArguments, const Kokkos::View<const double***>& p_fieldDofsCoordinates,
			 Kokkos::View<double***>& p_fieldDofsValues, const std::vector<int>& p_components,
			 const std::vector<const NSInputFunction::InputFunction*>& p_expressions ) {
  // Loop over the nodes of each entity and evaluate the prescribed condition on that node
  auto num_entity = p_bucket.GetSize();
  for( size_t e=0; e<num_entity; ++e ) {
    for( unsigned n=0; n<p_fieldNumNodesPerElement; ++n ) {
      // Copy the coordinates of that node
      for( int d=0; d<p_coordinatesNumDofsPerNode; ++d ) {
	p_functionArguments[d] = p_fieldDofsCoordinates( e, n, d );
      }
      
      // Evaluate the prescribed function
      if( p_components.size() == 0 ) {
	// This means that all the dofs of the node are prescribed
	// Get the pointer to the prescribed function
	auto input_function = p_expressions[0];
	
	// Loop over each dof defined on the node and evaluate the expression of the prescribed function
	for( int d=0; d<p_fieldNumDofsPerNode; ++d ) {
	  p_fieldDofsValues( e, n, d ) = input_function->Evaluate( p_functionArguments, d );
	}
      }
      else {
	// This means that some dofs of the node are prescribed
	// Loop over each specified component
	for( unsigned d=0; d<p_components.size(); ++d ) {
	  // Get the pointer to the input function
	  auto input_function = p_expressions[d];
	  
	  // Evaluate and set the value of the prescribed condition
	  // We set the second argument of Evaluate to 0 because the parsed function consists only of a single expression in this case
	  p_fieldDofsValues( e, n, p_components[d] ) = input_function->Evaluate( p_functionArguments, 0 );
	}
      }
    }
  }
}

      
void ApplyBoundaryConditionsOnFields( const NSMesh::Mesh* p_mesh, const NSField::Field<double>* p_coordinates,
				      const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
				      const std::map<std::string, const NSPrescribedCondition::PrescribedCondition*>& p_nameToBoundaryConditionMap,
				      double p_time ) {
  // Loop over each boundary condition
  std::vector<std::string> updated_fields;
  for( const auto& boundary_condition_info : p_nameToBoundaryConditionMap ) {
    // Get the pointer to the boundary condition
    const auto& boundary_condition = boundary_condition_info.second;

    // Get the expression of the boundary condition
    const auto& boundary_condition_expression = boundary_condition_info.second->GetPrescribedConditionFunction();

    // Get the vector of components of the field on which the boundary condition is specified
    const auto& components = boundary_condition_info.second->GetComponents();
    
    // Get the name of the field on which this condition is enforced
    const auto& field_name = boundary_condition->GetFieldName();
    updated_fields.push_back( field_name );
    
    // Get the name of the parts on which this condition is enforced
    const auto& boundary_parts_name = boundary_condition->GetPartsName();

    // Get the type of the closure entity (i.e. edge/face if 2D/3D)
    const auto& closure_entity_topology = p_mesh->GetClosureTopology();
    
    // Get the pointer to the field on which this condition is enforced
    NSExceptions::RuntimeError( p_nameToFieldMap.find( field_name ) == p_nameToFieldMap.end(), "The field on which a boundary condition is being enforced is not part of the finite element problem" );
    auto& field = p_nameToFieldMap.find( field_name )->second;

    // Get the number of dofs per node of the field
    auto field_num_dofs_per_node = field->GetNumDofsPerNode();

    // Get the data access type
    auto closure_access_type = ( p_mesh->SpatialDimension() == 2 ) ? NSField::AccessPer::ad_EDGE : NSField::AccessPer::ad_FACE;
    
    // Check that the boundary condition function is consistent with the number of dofs per node of the field
    if( components.size() > 0 ) {
      for( const auto& component : components ) {
	NSExceptions::InvalidArgument( component >= field_num_dofs_per_node, "A boundary condition is prescribed on a component that is not defined for the field" );
      }
    }
    else {
      NSExceptions::InvalidArgument( boundary_condition_expression[0]->GetNumComponents() != field_num_dofs_per_node, "The boundary condition function must have the same order as the field whose value it is prescribing" );
    }
    
    // Loop over the boundary mesh parts and enforce the boundary condition
    for( const auto& boundary_part_name : boundary_parts_name ) {
      // Get the pointer to the boundary mesh part
      auto boundary_mesh_part = p_mesh->GetBoundaryMeshPart( boundary_part_name );

      // Get the names of the interior mesh parts attached to this boundary mesh part
      const auto& attached_interior_mesh_parts_name = boundary_mesh_part->GetAttachedInteriorMeshPartsNames();
      
      // Loop over the interior mesh parts attached to the boundary mesh part
      for( const auto& interior_part_name : attached_interior_mesh_parts_name ) {
	// Get the buckets holding the closure entities of the boundary mesh part that are attached to the current interior mesh part
	std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> > closure_entity_buckets;
	p_mesh->GetBucketsByBoundaryPart( closure_entity_buckets, {NSMesh::Ownership::o_OWNED, NSMesh::Ownership::o_SHARED, NSMesh::Ownership::o_UNIVERSAL},
					  closure_entity_topology, boundary_part_name, interior_part_name );

	// Get the buckets holding the elements attached to the boundary mesh part and that belong  the current interior mesh part
	std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> > elements_buckets;
	p_mesh->GetBucketsByBoundaryPart( elements_buckets, {NSMesh::Ownership::o_OWNED, NSMesh::Ownership::o_SHARED, NSMesh::Ownership::o_UNIVERSAL},
					  NSMesh::Topology::tp_ELEMENT, boundary_part_name, interior_part_name );

	// Create the data structure that will hold the coordinates of the vertices of the closure entities of the current interior mesh part
	auto coordinate_field_element_type = p_coordinates->GetFunctionSpace( interior_part_name )->GetFieldElement();
	auto coordinates_num_nodes_per_closure_entity = NSField::factory_get_num_nodes_per_entity::Create( coordinate_field_element_type )( closure_entity_topology );
	auto coordinates_num_dofs_per_node = p_coordinates->GetNumDofsPerNode();
	Kokkos::View<double***> vertices_coordinates( "coordinates_boundary_conditions", BUCKETSIZE,
						      coordinates_num_nodes_per_closure_entity, coordinates_num_dofs_per_node );

	// Create the data structure that will hold the coordinates of the nodes defined on the closure entities attached to the current interior mesh part
	auto field_field_element_type = field->GetFunctionSpace( interior_part_name )->GetFieldElement();
	auto field_num_nodes_per_closure_entity = NSField::factory_get_num_nodes_per_entity::Create( field_field_element_type )( closure_entity_topology );
	Kokkos::View<double***> field_nodes_coordinates( "field_coordinates_boundary_conditions", BUCKETSIZE,
							 field_num_nodes_per_closure_entity, coordinates_num_dofs_per_node );
	
	// Loop over the buckets vector
	for( size_t b=0; b<closure_entity_buckets.size(); ++b ) {
	  // Get the closure entities bucket
	  const auto& closure_entity_bucket = closure_entity_buckets[b];

	  // Get the attached elements bucket
	  const auto& elements_bucket = elements_buckets[b];
	  
	  // Gather the coordinates of the corner nodes of the closure entitites attached to the current mesh part
	  p_coordinates->GatherData( closure_entity_bucket, p_coordinates->GetFunctionSpace( interior_part_name ), closure_access_type, vertices_coordinates );

	  // Get the coordinates of all the dofs of the field
	  NSField::factory_get_coordinates::Create( field_field_element_type )( closure_entity_bucket, p_mesh, closure_access_type,
										vertices_coordinates, field_nodes_coordinates );

	  // Create the vector to hold the arguments of the input function
	  // By default we assume that unknowns of the input function are x,y,z and t
	  std::vector<double> function_arguments = { 0.0, 0.0, 0.0, p_time };

	  // Create the Kokkos::View that will hold the evaluation of the boundary condition which will be later scattered to the field
	  Kokkos::View<double***> field_dofs_bc_values( "field_dofs_boundary_condition_values", BUCKETSIZE,
							field_num_nodes_per_closure_entity, field_num_dofs_per_node );
	  
	  // Set values per bucket
	  SetValuesPerBucket( elements_bucket, field_num_nodes_per_closure_entity, coordinates_num_dofs_per_node, field_num_dofs_per_node,
			      function_arguments, field_nodes_coordinates, field_dofs_bc_values, components, boundary_condition_expression );
 
	  // Scatter the field data after the boundary condition has been set
	  field->ScatterData( closure_entity_bucket, field->GetFunctionSpace( interior_part_name ), closure_access_type, field_dofs_bc_values );
	}
      }
    }
  }

  // Parallel communicate the shared data between the processors for each field that has been updated by a prescribed condition
  std::sort( updated_fields.begin(), updated_fields.end() );
  auto it = std::unique( updated_fields.begin(), updated_fields.end() );
  updated_fields.erase( it, updated_fields.end() );
  for( const auto& field_name : updated_fields ) {
    p_nameToFieldMap.find( field_name )->second->CommunicateInterfaceData( {NSMesh::Ownership::o_SHARED}, NSField::CombineType::c_REPLACE );
  }
}


void ApplyInitialConditionsOnFields( const NSMesh::Mesh* p_mesh, const NSField::Field<double>* p_coordinates,
				     const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
				     const std::map<std::string, const NSPrescribedCondition::PrescribedCondition*>& p_nameToInitialConditionMap,
				     double p_time ) {
  // Loop over each initial condition
  std::vector<std::string> updated_fields;
  for( const auto& initial_condition_info : p_nameToInitialConditionMap ) {
    // Get the pointer to the initial condition
    const auto& initial_condition = initial_condition_info.second;

    // Get the expression of the initial condition
    const auto& initial_condition_expression = initial_condition_info.second->GetPrescribedConditionFunction();

    // Get the vector of components of the field on which the initial condition is specified
    const auto& components = initial_condition_info.second->GetComponents();
    
    // Get the name of the field on which this condition is enforced
    const auto& field_name = initial_condition->GetFieldName();
    updated_fields.push_back( field_name );
    
    // Get the name of the parts on which this condition is prescribed
    const auto& interior_parts_name = initial_condition->GetPartsName();

    // Get the pointer to the field on which this condition is enforced
    NSExceptions::RuntimeError( p_nameToFieldMap.find( field_name ) == p_nameToFieldMap.end(), "The field on which a initial condition is being enforced is not part of the finite element problem" );
    auto& field = p_nameToFieldMap.find( field_name )->second;

    // Get the number of dofs per node
    auto field_num_dofs_per_node = field->GetNumDofsPerNode();
    
    // Check that the initial condition function is consistent with the number of dofs per node of the field
    if( components.size() > 0 ) {
      for( const auto& component : components ) {
	NSExceptions::InvalidArgument( component >= field_num_dofs_per_node, "An initial condition is prescribed on a component that is not defined for the field" );
      }
    }
    else {
      NSExceptions::InvalidArgument( initial_condition_expression[0]->GetNumComponents() != field_num_dofs_per_node, "The initial condition function must have the same order as the field whose values it is prescribing" );
    }
    
    // Loop over the interior mesh parts and enforce the initial condition
    for( const auto& interior_part_name : interior_parts_name ) {
      // Get the buckets holding the elements defined on the interior mesh part
      std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> > elements_buckets;
      p_mesh->GetBucketsByInteriorPart( elements_buckets, {NSMesh::Ownership::o_OWNED, NSMesh::Ownership::o_SHARED, NSMesh::Ownership::o_UNIVERSAL},
					NSMesh::Topology::tp_ELEMENT, interior_part_name );

      // Create the data structure that will hold the coordinates of the element vertices
      auto coordinate_field_element_type = p_coordinates->GetFunctionSpace( interior_part_name )->GetFieldElement();
      auto coordinates_num_nodes_per_element = NSField::factory_get_num_nodes_per_entity::Create( coordinate_field_element_type )( NSMesh::Topology::tp_ELEMENT );
      auto coordinates_num_dofs_per_node = p_coordinates->GetNumDofsPerNode();
      Kokkos::View<double***> vertices_coordinates( "coordinates_initial_conditions", BUCKETSIZE, coordinates_num_nodes_per_element, coordinates_num_dofs_per_node );

      // Create the data structure that will hold the coordinates of all the nodes of field
      auto field_field_element_type = field->GetFunctionSpace( interior_part_name )->GetFieldElement();
      auto field_num_nodes_per_element = NSField::factory_get_num_nodes_per_entity::Create( field_field_element_type )( NSMesh::Topology::tp_ELEMENT );
      Kokkos::View<double***> field_nodes_coordinates( "field_coordinates_initial_conditions", BUCKETSIZE,
						       field_num_nodes_per_element, coordinates_num_dofs_per_node );

      // Create the vector to hold the arguments of the input function
      // By default we assume that unknowns of the input function are x,y,z and t
      std::vector<double> function_arguments = { 0.0, 0.0, 0.0, p_time };
      
      // Loop over the buckets vector
      for( size_t b=0; b<elements_buckets.size(); ++b ) {
	// Get the attached elements bucket
	const auto& elements_bucket = elements_buckets[b];
	  
	// Gather the coordinates of the vertices of the element
	p_coordinates->GatherData( elements_bucket, p_coordinates->GetFunctionSpace( interior_part_name ), NSField::AccessPer::ad_ELEMENT, vertices_coordinates );

	// Get the coordinates of all the nodes of the field
	NSField::factory_get_coordinates::Create( field_field_element_type )( elements_bucket, p_mesh, NSField::AccessPer::ad_ELEMENT,
									      vertices_coordinates, field_nodes_coordinates );

	// Create the Kokkos::View that will hold the evaluation of the initial condition
	Kokkos::View<double***> field_dofs_values( "field_dofs_initial_condition_values", BUCKETSIZE,
						   field_num_nodes_per_element, field_num_dofs_per_node );

	// Evaluate the prescribed condition for the dofs defined at the nodes of the entities in the current bucket
	SetValuesPerBucket( elements_bucket, field_num_nodes_per_element, coordinates_num_dofs_per_node, field_num_dofs_per_node,
			    function_arguments, field_nodes_coordinates, field_dofs_values, components, initial_condition_expression );

	// Scatter the field data after the initial condition has been set
	field->ScatterData( elements_bucket, field->GetFunctionSpace( interior_part_name ), NSField::AccessPer::ad_ELEMENT, field_dofs_values );
      }
    }
  }

  // Parallel communicate the shared data between the processors for each field that has been updated by a prescribed condition
  std::sort( updated_fields.begin(), updated_fields.end() );
  auto it = std::unique( updated_fields.begin(), updated_fields.end() );
  updated_fields.erase( it, updated_fields.end() );
  for( const auto& field_name : updated_fields ) {
    p_nameToFieldMap.find( field_name )->second->CommunicateInterfaceData( {NSMesh::Ownership::o_SHARED}, NSField::CombineType::c_REPLACE );
  }
}


void ApplyBoundaryConditionsOnLinearSystem( const NSMesh::Mesh* p_mesh, const std::map<std::string, std::unique_ptr<NSField::Field<int> > >& p_nameToFieldDofIdsMap,
					    const std::map<std::string, const NSPrescribedCondition::PrescribedCondition*>& p_nameToBoundaryConditionMap,
					    NSLinearSystem::LinearSystem* const p_linearSystem ) {  
  // Loop over each boundary condition
  for( const auto& boundary_condition_info : p_nameToBoundaryConditionMap ) {
    // Get the pointer to the boundary condition
    const auto& boundary_condition = boundary_condition_info.second;

    // Get the vector of components of the field on which the boundary condition is specified
    const auto& components = boundary_condition_info.second->GetComponents();
    
    // Get the name of the field on which this condition is enforced
    const auto& field_name = boundary_condition->GetFieldName();
    
    // Get the pointer to the field holding the ids of the dofs of the field on which the boundary condition is enforced
    const auto& field_dof_vector = p_nameToFieldDofIdsMap.find( field_name )->second;

    // Get the number of dofs per node for the field
    auto field_num_dofs_per_node = field_dof_vector->GetNumDofsPerNode();
    
    // Get the name of the parts on which this condition is enforced
    const auto& boundary_parts_name = boundary_condition->GetPartsName();

    // Get the data access type
    auto closure_access_type = ( p_mesh->SpatialDimension() == 2 ) ? NSField::AccessPer::ad_EDGE : NSField::AccessPer::ad_FACE;
    
    // Loop over the boundary mesh parts and enforce the boundary condition
    for( const auto& boundary_part_name : boundary_parts_name ) {
      // Get the pointer to the boundary mesh part
      auto boundary_mesh_part = p_mesh->GetBoundaryMeshPart( boundary_part_name );

      // Get the names of the interior mesh parts attached to this boundary mesh part
      const auto& attached_interior_mesh_parts_name = boundary_mesh_part->GetAttachedInteriorMeshPartsNames();
      
      // Loop over the interior mesh parts attached to the boundary mesh part
      for( const auto& interior_part_name : attached_interior_mesh_parts_name ) {
	// Get the buckets holding the closure entities of the boundary mesh part that are attached to the current interior mesh part
	std::vector<NSMesh::NSEntity::Bucket<const NSMesh::NSEntity::Entity*> > closure_entity_buckets;
	p_mesh->GetBucketsByBoundaryPart( closure_entity_buckets, {NSMesh::Ownership::o_OWNED, NSMesh::Ownership::o_UNIVERSAL},
					  p_mesh->GetClosureTopology(), boundary_part_name, interior_part_name );

	// Loop over the buckets vector
	for( size_t b=0; b<closure_entity_buckets.size(); ++b ) {
	  // Get the closure entities bucket
	  const auto& closure_entity_bucket = closure_entity_buckets[b];

	  // Gather the dofs on which the prescribed conditions is enforced
	  auto field_field_element_type = field_dof_vector->GetFunctionSpace( interior_part_name )->GetFieldElement();
	  auto field_num_nodes_per_element = NSField::factory_get_num_nodes_per_entity::Create( field_field_element_type )( p_mesh->GetClosureTopology() );
	  Kokkos::View<int***> field_dofs_ids( "field_dofs_ids_boundary_condition", BUCKETSIZE,
					       field_num_nodes_per_element, field_num_dofs_per_node );
	  field_dof_vector->GatherData( closure_entity_bucket, field_dof_vector->GetFunctionSpace( interior_part_name ), closure_access_type, field_dofs_ids );
	  
	  // Create a vector holding the components of the field on which the boundary condition is applied
	  std::vector<int> field_bc_components = components;
	  unsigned dofs_offset = components.size();
	  if( components.size() == 0 ) {
	    field_bc_components.clear();
	    for( int d=0; d<field_num_dofs_per_node; ++d ) {
	      field_bc_components.push_back( d );
	    }
	    dofs_offset = field_num_dofs_per_node;
	  }
	  
	  // Loop over the closures entities in the bucket
	  auto num_entities = closure_entity_bucket.GetSize();
	  std::vector<int> all_dof_ids( num_entities * field_num_nodes_per_element * field_bc_components.size(), 0 ); // This vector is built for all entities in the bucket to be used for the rhs
	  std::vector<int> dof_ids_per_entity( field_num_nodes_per_element * field_bc_components.size(), 0 ); // This vector is built per entity to be used for the lhs
	  unsigned dofs_cnt = 0;	  
	  for( size_t e=0; e<num_entities; ++e ) {
	    for( unsigned n=0; n<field_num_nodes_per_element; ++n ) {
	      for( unsigned d=0; d<field_bc_components.size(); ++d ) {
		dof_ids_per_entity[ n * dofs_offset + d ] = field_dofs_ids( e, n, field_bc_components[d] );
		all_dof_ids[ dofs_cnt ] = field_dofs_ids( e, n, field_bc_components[d] );
		dofs_cnt++;
	      }
	    }
	    
	    // Set the diagonal of the lhs to 1 and the rest of the entries to zero at the rows corresponding to the current entity dofs
	    p_linearSystem->Diagonalize( dof_ids_per_entity.size(), dof_ids_per_entity.data(), 1.0 );
	  }
	  
	  // Set the rhs of the linear system to zero	  
	  p_linearSystem->SetRHSRowsToValue( all_dof_ids.size(), all_dof_ids.data(), 0.0 );
	}
      }
    }
  }
}
      
} // namespace NSFE
} // namespace NSProblem
} // namespace NSSimulation
