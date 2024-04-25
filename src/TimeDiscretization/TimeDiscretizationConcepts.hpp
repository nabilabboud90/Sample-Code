#ifndef TIMEDISCRETIZATIONCONCEPTS_HPP_
#define TIMEDISCRETIZATIONCONCEPTS_HPP_

#include <type_traits>

namespace NSTimeDiscretization {

// Define some constraints on the template parameters using C++ concepts
// Define general constraints on both the matrix type and the vector type
// Constraint to force the type to be a pointer
template<typename T>
concept Pointer = std::is_pointer_v<T>;

// Constraint to force the type to overload operator+=
template<typename T>
concept PlusEqual = requires ( T p_a, T p_b ) {
  *p_a += p_b;
};

// Constraint to force the type to overload operator-=
template<typename T>
concept MinusEqual = requires ( T p_a, T p_b ) {
  *p_a -= p_b;
};

// Constraint to force the type to overload operator*=
template<typename T>
concept TimesEqualDouble = requires ( T p_a, double p_b ) {
  *p_a *= p_b;
};

// Define the constraint on the vector type
template<typename VectorType>
concept VectorConstraint = Pointer<VectorType> && PlusEqual<VectorType> && MinusEqual<VectorType> && TimesEqualDouble<VectorType>;

// Define the constraint on the matrix type
template<typename MatrixType>
concept MatrixConstraint = Pointer<MatrixType> && PlusEqual<MatrixType> && MinusEqual<MatrixType> && TimesEqualDouble<MatrixType>;
  
// Define constraints that are specific to the matrix type
template<typename MatrixType, typename VectorType>
concept MatrixVectorMultiplication = requires ( MatrixType p_matrix, VectorType p_vector ) {
  *p_matrix * p_vector;
};

template<typename MatrixType, typename VectorType>
concept TimeDiscretizationParamterConstraints = MatrixVectorMultiplication<MatrixType, VectorType> && MatrixConstraint<MatrixType> && VectorConstraint<VectorType>;
  
} // namespace NSTimeDiscretization

#endif
