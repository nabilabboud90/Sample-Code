#include "NonlinearSolver/NonlinearSolvers.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "Factory/Factory.hpp"
#include "NonlinearSolver/NewtonSolver.hpp"

namespace NSNonlinearSolver {

NonlinearSolvers::NonlinearSolvers() {
}


void
NonlinearSolvers::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  p_datafile.Get( p_sectionPath + "/" + p_sectionName + "/LIST", _nonlinearSolversName );

  for( size_t ns=0; ns<_nonlinearSolversName.size(); ++ns ) {
    // Get the type of the nonlinear solver
    std::string nonlinear_solver_type = p_datafile.Get<std::string>( p_sectionPath + "/" + p_sectionName + "/" + _nonlinearSolversName[ns] + "/TYPE" );

    // Create the nonlinear solver
    std::unique_ptr<NSNonlinearSolver::NonlinearSolver> nonlinear_solver;
    nonlinear_solver.reset( NSNonlinearSolver::NonlinearSolver::factory_nonlinearsolver_Type::Create( nonlinear_solver_type ) );
    nonlinear_solver->Setup( p_datafile, p_sectionPath + "/" + p_sectionName, _nonlinearSolversName[ ns ] );

    // Store the nonlinear solver
    _nameToNonlinearSolverMap.insert( std::make_pair( _nonlinearSolversName[ ns ], std::move( nonlinear_solver ) ) );
  }
}


void
NonlinearSolvers::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  for( auto& nonlinear_solver : _nameToNonlinearSolverMap ) {
    nonlinear_solver.second->Initialize( p_com );
  }
}


const NSNonlinearSolver::NonlinearSolver*
NonlinearSolvers::GetNonlinearSolver( const std::string& p_name ) const {
  NSExceptions::InvalidArgument( _nameToNonlinearSolverMap.find( p_name ) == _nameToNonlinearSolverMap.end(), "The nonlinear solver you are asking for is not defined" );
  return _nameToNonlinearSolverMap.find( p_name )->second.get();
}

} // namespace LinearSolver
