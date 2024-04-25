#include "FEAlgorithm.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "Field/Fields.hpp"
#include "InputFunction/Functions.hpp"
#include "Mesh/Mesh.hpp"
#include "Equation/Equations.hpp"
#include "FELinearSystemHelper.hpp"
#include "LinearSystem/Epetra/EpetraLinearSystem.hpp"
#include "LinearSystem/LinearSystems.hpp"
#include "LinearSystem/LinearSolvers.hpp"
#include "FEPrescribedConditionsHelper.hpp"
#include "FESolutionToFieldDataHelper.hpp"
#include "Communicator/MPI/MPIComWrapper.hpp"
#include "NonlinearSolver/NonlinearSolver.hpp"
#include "Simulation/Problem/PDEProblem/DiscretizeInTime/DiscretizeInTime.hpp"

namespace NSSimulation { namespace NSProblem {

PDEAlgorithm* createFEAlgorithm() {
  return new FEAlgorithm;
}
    
FEAlgorithm::FEAlgorithm() {
}


void
FEAlgorithm::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  PDEAlgorithm::Setup( p_datafile, p_sectionPath, p_sectionName );

  // Parse the weak form information and create the weak form object
  // Read the names of the unknowns, test functions and input to the equation system
  p_datafile.Get( p_sectionPath + "/" + p_sectionName + "/EQUATIONSYSTEM/UNKNOWNS", _unknownsName );  
  p_datafile.Get( p_sectionPath + "/" + p_sectionName + "/EQUATIONSYSTEM/TESTFUNCTIONS", _testFunctionsName );
  p_datafile.Get( p_sectionPath + "/" + p_sectionName + "/EQUATIONSYSTEM/INPUTS", _inputsName );
  NSExceptions::InvalidArgument( _unknownsName.size() != _testFunctionsName.size(), "The number of UNKNOWNS and TESTFUNCTIONS in the input file must be the same" );
}


void
FEAlgorithm::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com,
			 const NSMesh::Mesh* const p_mesh, const NSField::Fields* const p_fields,
			 const NSPrescribedCondition::PrescribedConditions* const p_prescribedConditions,
			 const NSNonlinearSolver::NonlinearSolvers* const p_nonlinearSolvers,
			 const NSInputFunction::Functions* const p_functions,
			 const NSEquation::Equations* const p_equations,
			 const NSLinearSystem::LinearSystems* const p_linearSystems,
			 const NSLinearSystem::LinearSolvers* const p_linearSolvers ) {

  PDEAlgorithm::Initialize( p_com, p_mesh, p_fields, p_prescribedConditions, p_nonlinearSolvers, p_functions, p_equations, p_linearSystems, p_linearSolvers );
  
  // Set problem unknown and test functions
  _initializeProblemFields( p_fields );
  
  // Set weak form
  _initializeWeakForm( p_equations, p_mesh->SpatialDimension() );
  
  // Set the linear system
  _initializeLinearSystem( p_com, p_linearSystems, p_linearSolvers, p_fields );

  // Apply the initial conditions on fields if any
  NSFE::ApplyInitialConditionsOnFields( _mesh, _nameToFieldMap.find( "COORDINATES" )->second, _nameToFieldMap, _nameToICMap, 0.0 );

  // Apply the boundary conditions on the fields if any
  NSFE::ApplyBoundaryConditionsOnFields( _mesh, _nameToFieldMap.find( "COORDINATES" )->second, _nameToFieldMap, _nameToStrongBCMap, 0.0 );

  // Initialize the solution vector
  NSFE::SetSolutionVectorFromFieldsData( _mesh, _nameToFieldMap, _nameToFieldDofIdsMap, _linearSystem.get() );

  // Initialize the finite element assembly operator whose in charge of assembling the lhs/rhs
  _assemblyOperator = std::make_unique<NSFE::FEAssemblyOperator>( true );
  _assemblyOperator->Initialize( _mesh, _unknownsName, _testFunctionsName, _inputsName, _weakForms.get(),
				 _nameToFieldMap, _nameToIdMap, _nameToFieldDofIdsMap, _linearSystem->GetLhs() );
}


void
FEAlgorithm::_initializeProblemFields( const NSField::Fields* const p_fields ) {
  // Set the pointers to the unknowns ofthe problem
  unsigned test_function_cnt = 0;
  unsigned all_fields_cnt = 0;
  for( const auto& unknown_name : _unknownsName ) {
    _nameToFieldMap.insert( std::make_pair( unknown_name, &p_fields->GetField<double>( unknown_name ) ) );
    _nameToIdMap.insert( std::make_pair( unknown_name, test_function_cnt ) );
    all_fields_cnt++;
    
    _testFunctionNameToOrderMap.insert( std::make_pair( _testFunctionsName[ test_function_cnt ], p_fields->GetField<double>( unknown_name ).GetOrder() ) );
    _nameToIdMap.insert( std::make_pair( _testFunctionsName[ test_function_cnt ], _unknownsName.size() + test_function_cnt ) );
    all_fields_cnt++;
    
    test_function_cnt++;
  }
  
  // Add the input fields
  for( const auto& input_name : _inputsName ) {    
    _nameToFieldMap.insert( std::make_pair( input_name, &p_fields->GetField<double>( input_name ) ) );
    _nameToIdMap.insert( std::make_pair( input_name, all_fields_cnt ) );
    all_fields_cnt++;
  }
  
  // Add the coordinates field
  _nameToFieldMap.insert( std::make_pair( "COORDINATES", &p_fields->GetField<double>( "COORDINATES" ) ) );
  _nameToIdMap.insert( std::make_pair( "COORDINATES", all_fields_cnt ) );
}

    
void
FEAlgorithm::_initializeWeakForm( const NSEquation::Equations* const p_equations, int p_ndim ) {
  // Create and intialize the weak form object
  const auto& equation_system = p_equations->GetEquationSystem( _equationSystemName );
  _weakForms = std::make_unique<NSEquation::NSWeakForm::WeakForms>( equation_system->GetEquationNames(), equation_system->GetEquationExpressions() );
  
  // Create a map from the name of the fields involved in the weak form to the order of the fields
  std::unordered_map<std::string, std::pair<int, NSTypes::FieldOrder> > field_to_order_map;

  // Add the unknowns
  for( const auto& name : _unknownsName ) {
    int field_id = _nameToIdMap.find( name )->second;
    field_to_order_map.insert( std::make_pair( name, std::make_pair( field_id, _nameToFieldMap.find( name )->second->GetOrder() ) ) );
  }

  // Add the test functions
  for( const auto& name : _testFunctionsName ) {
    int field_id = _nameToIdMap.find( name )->second;
    field_to_order_map.insert( std::make_pair( name, std::make_pair( field_id, _testFunctionNameToOrderMap.find( name )->second ) ) );
  }
  
  // Add the input fields order
  for( const auto& input_name : _inputsName ) {
    int field_id = _nameToIdMap.find( input_name )->second;
    field_to_order_map.insert( std::make_pair( input_name, std::make_pair( field_id, _nameToFieldMap.find( input_name )->second->GetOrder() ) ) );
  }

  // Initialize the weak form
  _weakForms->Initialize( field_to_order_map, p_ndim );
}


void
FEAlgorithm::_createLinearSystem( const NSLinearSystem::LinearSystems* const p_linearSystems, const NSLinearSystem::LinearSolvers* const p_linearSolvers ) {
  // Get the linear system data
  const auto& linear_system_data = p_linearSystems->GetLinearSystemData( _linearSystemName );
  
  // Get the name of the linear solver to solve this linear system
  const auto& linear_solver_name = linear_system_data->GetLinearSolverName();

  // Get the object holding data about the specified linear solver
  const auto& linear_solver_data = p_linearSolvers->GetLinearSolver( linear_solver_name );
  
  // Create the linear system
  _linearSystem.reset( NSLinearSystem::LinearSystem::factory_linearsystem_Type::Create( linear_system_data->GetType() ) );
  _linearSystem->Setup( linear_system_data.get(), linear_solver_data.get() );
}

    
void
FEAlgorithm::_initializeLinearSystem( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com,
				      const NSLinearSystem::LinearSystems* const p_linearSystems, const NSLinearSystem::LinearSolvers* const p_linearSolvers,
				      const NSField::Fields* const p_fields ) {

  // Create a map that holds the names of the mesh parts for which degrees of freedom of a given field will be generated
  std::map<int, std::vector<std::string> > field_id_to_mesh_parts_name;
  std::vector<int> fields_id;
  for( unsigned i=0; i<_unknownsName.size(); ++i ) {
    field_id_to_mesh_parts_name.insert( std::make_pair( i, std::vector<std::string>() ) );
    fields_id.push_back( i );
  }
  _weakForms->GetPartNamesPerUnknown( fields_id, NSEquation::NSWeakForm::WeakFormTermType::wf_INTERIOR, field_id_to_mesh_parts_name );

  // Set the degrees of freedom ids
  std::map<std::string, std::vector<std::string> > mesh_part_name_to_field_name;
  auto linear_system_block_structure = NSSimulation::NSProblem::NSFE::SetFEDofsIds( p_com, p_fields, mesh_part_name_to_field_name, _unknownsName,
										    field_id_to_mesh_parts_name, _nameToFieldMap, _nameToFieldDofsIdFs,
										    _nameToFieldDofIdsMap, _mesh, _linearSystemName );
  
  // Create the linear system object and parse its data
  _createLinearSystem( p_linearSystems, p_linearSolvers );
  
  // Build the vector holding all ids of the dofs
  unsigned total_num_owned_dofs = 0;
  std::vector<int> all_dofs_ids;
  for( const auto& unknown_name : _unknownsName ) {
    const auto& field_dofs = _nameToFieldDofIdsMap.find( unknown_name )->second;
    const auto& dof_ids = field_dofs->GetFieldData();
    all_dofs_ids.insert( all_dofs_ids.end(), dof_ids.begin(), dof_ids.end() );

    // Increment the total number of dofs defined on owned entities of this processor
    const auto& all_function_spaces = field_dofs->GetAllFunctionSpaces();
    for( const auto& function_space : all_function_spaces ) {
      total_num_owned_dofs += field_dofs->GetNumDofsPerNode() * function_space->GetTotalNumDofsPerOwnerhsip( {NSMesh::Ownership::o_OWNED} );
    }
  }
  
  // Initialize the linear system
  _linearSystem->Initialize( _com, all_dofs_ids, total_num_owned_dofs, NSFE::GetFEMatrixGraphBuilder( _mesh, mesh_part_name_to_field_name, _nameToFieldDofIdsMap ),
			     linear_system_block_structure );
}


void
FEAlgorithm::Execute() {
  auto linearize = [&]() {
    // Zero the lhs, rhs and the solution increment vector of the linear system
    _linearSystem->ZeroSystem();
    
    // Update the fields before the assembly of the linear system
    // This call is responsible for enforcing any strongly applied boundary conditions on the fields of the problem at the current step
    NSFE::ApplyBoundaryConditionsOnFields( _mesh, _nameToFieldMap.find( "COORDINATES" )->second, _nameToFieldMap, _nameToStrongBCMap, 0.0 );
    
    // Discretize in space
    _assemblyOperator->Assemble( _mesh, _nameToFieldMap, _nameToIdMap, _nameToFieldDofIdsMap, _linearSystem->GetVector( 0 ), _linearSystem->GetLhs() );
    
    // Parallel communicate the processors shared data of the global rhs 
    auto rhs = _linearSystem->GetVector( 0 );    
    rhs->GlobalAssemble( NSLinearSystem::CombineMode::cm_ADD );
    
    // Parallel communicate the processors shared data of the global lhs 
    auto lhs = _linearSystem->GetLhs();
    lhs->GlobalAssemble( NSLinearSystem::CombineMode::cm_ADD );
    
    // Discretize in time
    std::vector<double> dts( 1, 0.0 );
    DiscretizeInTime( _linearSystem.get(), _timeDiscretization, dts );
  };

  auto finalize = [&]() -> std::vector<double> {
    // Apply the strongly enforced boundary conditions to the lhs and rhs of the linear system
    NSFE::ApplyBoundaryConditionsOnLinearSystem( _mesh, _nameToFieldDofIdsMap, _nameToStrongBCMap, _linearSystem.get() );
    
    // Get the inf norma and the scaled l2 norm of the rhs
    auto rhs = _linearSystem->GetVector( 0 );
    std::vector<double> norm( 2, 0.0 );
    norm[0] = rhs->NormInf();
    norm[1] = rhs->Norm2() / rhs->GetGlobalSize();
    
    return norm;
  };

  auto solve = [&]() {
    // Solve the linear system
    _linearSystem->Solve();
    
    // Update the unknown fields with the solution increment
    NSFE::UpdateFieldDataFromSolutionIncrement( _mesh, _nameToFieldMap, _nameToFieldDofIdsMap, _linearSystem.get() );
  };

  auto print = [&]( double p_scaledL2Norm, double p_infinityNorm, int p_iterationNumber, double p_tolerance ) {
    if( _com->GetRank() == 0 ) {
      std::cout << std::endl;
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cout << "      Nonlinear solver:      linear system          =      " << _linearSystemName << std::endl;
      std::cout << "                             iteration              =      " << p_iterationNumber << std::endl;
      std::cout << "                             scaledResidualL2Norm   =      " << p_scaledL2Norm << std::endl;
      std::cout << "                             residualInfNorm        =      " << p_infinityNorm << std::endl;
      std::cout << "                             tolerance              =      " << p_tolerance << std::endl;
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    }
  };
  
  _nonlinearSolver->Iterate( linearize, finalize, solve, print, true );
}
    
} // namespace NSProblem
} // namespace NSSimulation
