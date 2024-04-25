#include "DiscretizeInTime.hpp"
#include "TimeDiscretization/TimeDiscretization.hpp"
#include "TimeDiscretization/Quasistatic.hpp"
#include "LinearSystem/Matrix.hpp"
#include "LinearSystem/Vector.hpp"
#include "LinearSystem/LinearSystem.hpp"

namespace NSSimulation { namespace NSProblem {

void DiscretizeInTimeQuasistatic( NSLinearSystem::Vector* const p_solution, NSLinearSystem::Matrix* const p_matrix,
				  NSLinearSystem::Vector* const p_residual, const std::vector<double>& p_dts ) {
  NSTimeDiscretization::TimeDiscretization<NSLinearSystem::Matrix*, NSLinearSystem::Vector*, NSTimeDiscretization::Scheme::ts_QUASISTATIC> time_integrator;

  // Build the vector of matrices and rhs/solutions to be passed to the quasistatic time integrator, if needed
  std::vector<NSLinearSystem::Matrix*> matrices;
  matrices.push_back( p_matrix );
    
  std::vector<NSLinearSystem::Vector*> vectors;
  vectors.push_back( p_residual );
  vectors.push_back( p_solution );
    
  // Complete the linear system
  time_integrator( matrices, vectors, p_dts );
}

    
void DiscretizeInTime( NSLinearSystem::LinearSystem* const p_linearSystem, const NSTimeDiscretization::Scheme& p_scheme, const std::vector<double>& p_dts ) {
  // Get the solution
  auto solution = p_linearSystem->GetVector( 1 );

  // Get the lhs
  auto matrix = p_linearSystem->GetLhs();

  // Get the rhs
  auto residual = p_linearSystem->GetVector( 0 );
  
  // Discretize in time
  if( p_scheme == NSTimeDiscretization::Scheme::ts_QUASISTATIC ) {
    DiscretizeInTimeQuasistatic( solution, matrix, residual, p_dts );
  }
}

} // namespace NSProblem
} // namespace NSSimulation
