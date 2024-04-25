#ifndef TIMEDISCRETIZATION_HPP_
#define TIMEDISCRETIZATION_HPP_

#include <vector>

#include "TimeDiscretizationEnum.hpp"
#include "Utility/Exceptions.hpp"
#include "TimeDiscretizationConcepts.hpp"

namespace NSTimeDiscretization {
  
template<typename MatrixType, typename VectorType, NSTimeDiscretization::Scheme Scheme> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
class TimeDiscretization {
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
   * For example, for a BDF1, this function must do the following
   * RHS = M * (Y_{n+1} - Y_{n}) / dt - K * Y_{n+1} - F                                                                                                                   
   * LHS = M * del{Y}_{n+1} / dt - K * del{Y}_{n+1}
   * M is the mass matrix, K is the stiffness matrix, Y is the solution vector, del{Y} is the increment in the solution vector
   @param matrices - A vector of all the matrices involved in the system, mainly the mass matrix and the stiffness matrix
   @param vectors - A vector holding all the vectors involved in the system, mainly the residual and the solution vectors at different discrete times
   @param dts - A vector holding the last N time step sizes of the simulation
  */
  void operator()( std::vector<MatrixType>& p_matrices, std::vector<VectorType>& p_vectors, const std::vector<double>& p_dts );
};

template<typename MatrixType, typename VectorType, NSTimeDiscretization::Scheme Scheme> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
TimeDiscretization<MatrixType, VectorType, Scheme>::TimeDiscretization() {
}


template<typename MatrixType, typename VectorType, NSTimeDiscretization::Scheme Scheme> requires TimeDiscretizationParamterConstraints<MatrixType, VectorType>
void
TimeDiscretization<MatrixType, VectorType, Scheme>::operator()( [[maybe_unused]]std::vector<MatrixType>& p_matrices, [[maybe_unused]]std::vector<VectorType>& p_vectors, [[maybe_unused]]const std::vector<double>& p_dts ) {
  NSExceptions::InvalidArgument( true, "The time discretization is not implemented for the template parameters used" );
}

} // namespace NSTimeDiscretization

#endif
