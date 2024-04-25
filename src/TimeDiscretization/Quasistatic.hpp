#ifndef QUASISTATIC_HPP_
#define QUASISTATIC_HPP_

#include <vector>

#include "TimeDiscretizationEnum.hpp"
#include "Utility/Exceptions.hpp"
#include "TimeDiscretizationConcepts.hpp"

namespace NSTimeDiscretization {

template<typename MatrixType, typename VectorType> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
class TimeDiscretization<MatrixType, VectorType, NSTimeDiscretization::Scheme::ts_QUASISTATIC> {
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
   * For Quasistatic time discretization, this function must do the following
   * RHS = K * Y_{n+1} - F, if the problem is linear
   * If the problem is nonlinear, then the RHS will be constructed by building element vectors and assembling them directly into the RHS
   @param matrices - A vector of all the matrices involved in the system, mainly the mass matrix and the stiffness matrix
   @param vectors - A vector holding all the vectors involved in the system, mainly the residual and the solution vectors at different discrete times
   @param dts - A vector holding the last N time step sizes of the simulation
  */
  void operator()( std::vector<MatrixType>& p_matrices, std::vector<VectorType>& p_vectors, const std::vector<double>& p_dts );
};

template<typename MatrixType, typename VectorType> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
TimeDiscretization<MatrixType, VectorType, NSTimeDiscretization::Scheme::ts_QUASISTATIC>::TimeDiscretization() {
}


template<typename MatrixType, typename VectorType> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
void
TimeDiscretization<MatrixType, VectorType, NSTimeDiscretization::Scheme::ts_QUASISTATIC>::operator()( std::vector<MatrixType>& p_matrices, std::vector<VectorType>& p_vectors, [[maybe_unused]]const std::vector<double>& p_dts ) {
  // The p_vectors must have two vectors in the following order
  // p_vector[0] must have a pointer to the residual vector
  // p_vector[1] must have a pointer to the solution vector
  NSExceptions::InvalidArgument( p_vectors.size() != 2, "The p_vectors passed to the Complete method of the TimeDiscretization partially specialized with ts_QUASISTATIC must be of size 2" );

  // The p_matrices must have a single matrix
  // The matrix, if needed, should be defined by the linear terms of the partial differential equation
  NSExceptions::InvalidArgument( p_matrices.size() != 1, "The p_matrices passed to the Complete method of the TimeDiscretization partially specialized with ts_QUASISTATIC must be of size 1" );

  // Get the residual vector
  auto residual = p_vectors[0];

  // Get the solution vector
  auto solution = p_vectors[1];

  // Get the matrix
  auto matrix = p_matrices[0];
  
  // Update the residual vector
  *residual += ((*matrix) * solution).get();
}

} // namespace NSTimeDiscretization

#endif
