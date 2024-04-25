#ifndef BDF1_HPP_
#define BDF1_HPP_

#include <vector>

#include "TimeDiscretizationEnum.hpp"
#include "Utility/Exceptions.hpp"
#include "TimeDiscretizationConcepts.hpp"

namespace NSTimeDiscretization {

template<typename MatrixType, typename VectorType> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
class TimeDiscretization<MatrixType, VectorType, NSTimeDiscretization::Scheme::ts_BDF1> {
public:
  /*
   * Default constructor
   */
  explicit TimeDiscretization();
  
  /*
   * Copy constructor not implemented
   @param ti - time discretization to copy from
  */
  TimeDiscretization( const TimeDiscretization& p_ti ) = delete;
  
  /*
   * Equal operator not implemented
   @param ti - time discretization to set equal to
  */ 
  TimeDiscretization& operator=( const TimeDiscretization& p_ti ) = delete;
  
  /*
   * Function that completes the lhs and rhs of a time-dependent ordinary differential equation
   * For BDF1 time discretization, this function must do the following
   // M * (Y_{n+1} - Y_{n}) / dt = K * Y_{n+1} + F
   // Linearizing => M * (Y_{n+1} - Y_{n}) / dt - K * Y_{n+1} - F + M * del{Y}_{n+1} / dt - K * del{Y}_{n+1} = 0
   // Define, RHS = M * (Y_{n+1} - Y_{n}) / dt - K * Y_{n+1} - F
   // and,    LHS = M * del{Y}_{n+1} / dt - K * del{Y}_{n+1}
   * If the problem is nonlinear, then the RHS will be constructed by building element vectors and assembling them directly into the RHS
   @param matrices - A vector of all the matrices involved in the system
   @param vectors - A vector holding all the vectors involved in the system
   @param dts - A vector holding the last N time step sizes of the simulation
  */
  void operator()( std::vector<MatrixType>& p_matrices, std::vector<VectorType>& p_vectors, const std::vector<double>& p_dts );
};

template<typename MatrixType, typename VectorType> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
TimeDiscretization<MatrixType, VectorType, NSTimeDiscretization::Scheme::ts_BDF1>::TimeDiscretization() {
}


template<typename MatrixType, typename VectorType> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
void
TimeDiscretization<MatrixType, VectorType, NSTimeDiscretization::Scheme::ts_BDF1>::operator()( std::vector<MatrixType>& p_matrices, std::vector<VectorType>& p_vectors, const std::vector<double>& p_dts ) {
  // The p_vectors must have three vectors in the following order
  // p_vector[0] must have a pointer to the residual vector
  // p_vector[1] must have a pointer to the solution vector at tnp1
  // p_vector[2] must have a pointer to the solution vector at tn
  NSExceptions::InvalidArgument( p_vectors.size() != 3, "The p_vectors passed to the Complete method of the TimeDiscretization partially specialized with ts_QUASISTATIC must be of size 3" );

  // The p_matrices must have at least two matrices
  // If the problem has both linear and nonlinear terms, then p_matrices must have three matrices
  // The matrices must be in the following order
  // p_matrices[0] must have a pointer to the system matrix
  // p_matrices[1] must have a pointer to the mass matrix
  // p_matrices[2] if needed, must have a pointer to the matrix holding the linear terms of the partial differential equation
  NSExceptions::InvalidArgument( p_matrices.size() < 2 || p_matrices.size() > 3, "The p_matrices passed to the Complete method of the TimeDiscretization partially specialized with ts_BDF1 must be of size 2 or size 3" );

  // Get the residual vector
  auto residual = p_vectors[0];

  // Get the solution vector
  const auto* const solution_np1 = p_vectors[1];
  const auto* const solution_n = p_vectors[2];
  
  // Get the system matrix
  auto system_matrix = p_matrices[0];
  
  // Get the mass matrix
  const auto* const mass_matrix = p_matrices[1];
  
  // The starting point is the following
  // M * (Y_{n+1} - Y_{n}) / dt = K * Y_{n+1} + F
  // Linearizing => M * (Y_{n+1} - Y_{n}) / dt - K * Y_{n+1} - F + M * del{Y}_{n+1} / dt - K * del{Y}_{n+1} = 0
  // Define, RHS = M * (Y_{n+1} - Y_{n}) / dt - K * Y_{n+1} - F
  // and,    LHS = M * del{Y}_{n+1} / dt - K * del{Y}_{n+1}
  // Implementation of the above algebraic system is below

  // Compute RHS
  // RHS = -F
  *residual *= -1.0;

  // RHS = -K * Y_{n+1} - F, where K is the linearization of the linear terms
  if( p_matrices.size() == 3 ) {
    const auto klMat = p_matrices[2];
    *residual -= *( std::move( (*klMat) * (*solution_np1) ).get() );
  }

  // RHS = ( -K * Y_{n+1} - F ) * dt
  *residual *= p_dts.back();

  // RHS = ( -K * Y_{n+1} - F ) * dt - M * Y_{n}
  *residual -= *( std::move( (*mass_matrix) * (*solution_n) ).get() );

  // RHS = ( -K * Y_{n+1} - F ) * dt - M * Y_{n} + M * Y_{n+1}
  *residual += *( std::move( (*mass_matrix) * (*solution_np1) ).get() );

  // Compute LHS
  // LHS = -K * dt
  *system_matrix *= ( -1.0 * p_dts.back() );

  // LHS = M - K * dt
  *system_matrix += *mass_matrix;

}

} // namespace NSTimeDiscretization

#endif
