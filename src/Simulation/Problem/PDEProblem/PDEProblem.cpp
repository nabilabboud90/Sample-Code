#include "PDEProblem.hpp"
#include "Mesh/STKMesh/STKMesh.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "PDEAlgorithms/FiniteElement/FEAlgorithm.hpp"
#include "Factory/Factory.hpp"

namespace NSSimulation { namespace NSProblem {

Problem* createPDEProblem() {
  return new PDEProblem;
}
      
      
PDEProblem::PDEProblem() {
}


void
PDEProblem::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  Problem::Setup( p_datafile, p_sectionPath, p_sectionName );
  
  // Create and setup the mesh
  std::string mesh_type = p_datafile.Get<std::string>( p_sectionPath + "/MESH/TYPE" );
  _mesh.reset( NSMesh::Mesh::factory_mesh_Type::Create( mesh_type ) );
  _mesh->Setup( p_datafile, p_sectionPath, "MESH" );
  
  // Create and setup the fields
  _fields = std::make_unique<NSField::Fields>();
  _fields->Setup( p_datafile, p_sectionPath, "FIELDS" );
  
  // Create and setup the input functions
  _functions = std::make_unique<NSInputFunction::Functions>();
  _functions->Setup( p_datafile, p_sectionPath, "FUNCTIONS" );
  
  // Create and setup the prescribed conditions
  _prescribedConditions = std::make_unique<NSPrescribedCondition::PrescribedConditions>();
  _prescribedConditions->Setup( p_datafile, p_sectionPath, "PRESCRIBEDCONDITIONS" );
  
  // Create and setup the nonlinear solvers
  _nonlinearSolvers = std::make_unique<NSNonlinearSolver::NonlinearSolvers>();
  _nonlinearSolvers->Setup( p_datafile, p_sectionPath, "NONLINEARSOLVERS" );
  
  // Create and setup the output
  _outputs = std::make_unique<NSOutput::Outputs>();
  _outputs->Setup( p_datafile, p_sectionPath, "OUTPUTS" );
  
  // Create and setup the linear solvers
  _linearSolvers = std::make_unique<NSLinearSystem::LinearSolvers>();
  _linearSolvers->Setup( p_datafile, p_sectionPath, "LINEARSOLVERS" );
  
  // Create and setup the linear systems
  _linearSystems = std::make_unique<NSLinearSystem::LinearSystems>();
  _linearSystems->Setup( p_datafile, p_sectionPath, "LINEARSYSTEMS" );
  
  // Create and setup the equations
  _equations = std::make_unique<NSEquation::Equations>();
  _equations->Setup( p_datafile, p_sectionPath, "EQUATIONS" );
  
  // Create and setup the algorithms used to solve the problem
  p_datafile.Get( p_sectionPath + "/SOLUTION/ALGORITHMS", _algorithmNames );
  for( const auto& algorithm_name : _algorithmNames ) {
    // Read the type of the algorithm
    std::string algorithm_type = p_datafile.Get<std::string>( p_sectionPath + "/SOLUTION/" + algorithm_name + "/TYPE" );

    // Create the algorithm
    std::unique_ptr<NSSimulation::NSProblem::PDEAlgorithm> pde_algorithm;
    pde_algorithm.reset( NSSimulation::NSProblem::PDEAlgorithm::factory_pde_algorithm_Type::Create( algorithm_type ) );

    // Setup the algorithm
    pde_algorithm->Setup( p_datafile, p_sectionPath + "/SOLUTION", algorithm_name );

    // Store the algorithm
    _nameToAlgorithmMap.insert( std::make_pair( algorithm_name, std::move( pde_algorithm ) ) );
  }
}


void
PDEProblem::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  Problem::Initialize( p_com );

  // Initialize the mesh
  _mesh->Initialize( p_com );
  
  // Initialize the fields
  _fields->Initialize( p_com, _mesh.get() );
  
  // Initialize the input functions
  _functions->Initialize( _mesh->SpatialDimension() );
  
  // Initialize the prescribed conditions
  _prescribedConditions->Initialize( _functions.get(), _mesh->SpatialDimension() );
  
  // Initialize the nonlinear solvers
  _nonlinearSolvers->Initialize( p_com );
  
  // Initialize the outputs
  _outputs->Initialize( p_com, _mesh.get(), _fields.get() );
  
  // Initialize the pde algorithms
  for( const auto& algorithm_name : _algorithmNames ) {
    _nameToAlgorithmMap.find( algorithm_name )->second->Initialize( p_com, _mesh.get(), _fields.get(), _prescribedConditions.get(),
								    _nonlinearSolvers.get(), _functions.get(), _equations.get(),
								    _linearSystems.get(), _linearSolvers.get() );
  }

  _outputs->WriteOutput( 0.0 );
}


void
PDEProblem::Execute() {
  for( auto& pde_algorithm : _nameToAlgorithmMap ) {
    pde_algorithm.second->Execute();
  }

  _outputs->WriteOutput( 1.0 );
  _outputs->Finalize();
}
    
} // namespace NSProblem
} // namespace NSSimulation
