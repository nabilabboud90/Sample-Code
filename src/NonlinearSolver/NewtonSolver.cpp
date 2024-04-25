#include "NonlinearSolver/NewtonSolver.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "Communicator/MPI/MPIComWrapper.hpp"

namespace NSNonlinearSolver {

NonlinearSolver* createNewtonSolver() {
  return new NewtonSolver;
}
  
  
NewtonSolver::NewtonSolver() {
}

  
void
NewtonSolver::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  NonlinearSolver::Setup( p_datafile, p_sectionPath, p_sectionName );

  // Read the value of eta maximum
  _etaMaximum = p_datafile.Get<double>( p_sectionPath + "/" + p_sectionName + "/ETAMAXIMUM" );
}


void
NewtonSolver::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  NonlinearSolver::Initialize( p_com );
}


bool
NewtonSolver::Iterate( std::function<void()> p_linearize, std::function<std::vector<double>()> p_finalize, std::function<void()> p_solve, std::function<void( double, double, int, double )> p_print, bool p_verbose ) const {
  // Create some scalars that will be needed in the newton loop
  int iterationNumber = 0;
  double ratio = 0.0;
  std::vector<double> previousResidualNorms( 2, 1.0 );
  std::vector<double> residualNorms( 2, 1.0 );
  double scaledL2ResidualNorm = 0.0;
  double firstResidualL2Norm = 0.0;
  double linearRelTol = std::fabs( _etaMaximum );
  double eta_old = 0.0;
  double eta_new = 0.0;

  // Newton loop
  while( iterationNumber < _maxIterations ) {
    ratio = residualNorms[0] / previousResidualNorms[0];
    previousResidualNorms[0] = residualNorms[0];

    // Call the assemble and update of the rhs and lhs of the linearized problem
    p_linearize();
    
    // Finalize the linear system, i.e, do some parallel communication and compute some norms of the residual
    residualNorms = p_finalize();

    // Store the L2 norm of the residual at the 1st iteration to be used later for normalization
    if( iterationNumber == 0 ) {
      firstResidualL2Norm = residualNorms[1];
    }
    scaledL2ResidualNorm = residualNorms[1] / std::max( std::numeric_limits<double>::epsilon(), firstResidualL2Norm );

    // Print some output to screen about the newton iteration
    if( p_verbose ) {
      p_print( scaledL2ResidualNorm, residualNorms[0], iterationNumber, _absoluteTolerance );
    }

    // Check convergence
    if( ( ( scaledL2ResidualNorm < _absoluteTolerance ) || ( residualNorms[0] < _absoluteTolerance ) ) && iterationNumber > 0 ) {
      return true;
    }

    // If not converged solve the linear system
    p_solve();

    // Increment the iteration number
    iterationNumber++;

    // This is needed for the linear solver
    if( _etaMaximum > 0 ) {
      eta_old = linearRelTol;
      eta_new = _gamma * ratio * ratio;
      
      if( _gamma * eta_old * eta_old > 0.1 ) {
	eta_new = std::max<double>( eta_new, _gamma * eta_old * eta_old );
      }

      linearRelTol = std::min<double>( eta_new, _etaMaximum );
      linearRelTol = std::min<double>( _etaMaximum, std::max<double> ( linearRelTol, 0.5 * _absoluteTolerance / residualNorms[0] ) );
    }
  }

  return false;
}

} // namespace NonlinearSolver
