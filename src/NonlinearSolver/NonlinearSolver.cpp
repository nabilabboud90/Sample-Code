#include "NonlinearSolver/NonlinearSolver.hpp"
#include "Parser/InputParser/Datafile.hpp"

namespace NSNonlinearSolver {

NonlinearSolver::NonlinearSolver() {
}


void
NonlinearSolver::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  // Set the name of the nonlinear solver
  _name = p_sectionName;

  // Set the maximum number of iterations
  _maxIterations = p_datafile.Get<int>( p_sectionPath + "/" + p_sectionName + "/MAXITERATIONS" );

  // Set the relative tolerance
  _relativeTolerance = p_datafile.Get<double>( p_sectionPath + "/" + p_sectionName + "/RELATIVETOLERANCE" );

  // Set the absolute tolerance
  _absoluteTolerance = p_datafile.Get<double>( p_sectionPath + "/" + p_sectionName + "/ABSOLUTETOLERANCE" );
}
  
  
void
NonlinearSolver::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  _com = p_com;
}
  
} // namespace NSNonlinearSolver
