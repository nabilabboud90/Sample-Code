#ifndef FEALGORITHM_HPP_
#define FEALGORITHM_HPP_

#include "Simulation/Problem/PDEProblem/PDEAlgorithms/PDEAlgorithm.hpp"
#include "Factory/Factory.hpp"
#include "Equation/WeakForm/WeakForms.hpp"
#include "FEAssemblyOperator.hpp"

#include <Kokkos_Core.hpp>

namespace NSField { template<typename Type> class Field; }
namespace NSField { class FunctionSpace; }

namespace NSSimulation { namespace NSProblem {
    
class FEAlgorithm : public PDEAlgorithm {
public:
  /*
   * Default constructor
   */
  explicit FEAlgorithm();

  /*
   * Copy constructor not implemented
   @param feAlgorithm - FEAlgorithm object used to copy from
   */
  FEAlgorithm( const FEAlgorithm& p_feAlgorithm ) = delete;

  /*
   * Equal operator not implemented
   @param feAlgorithm - FEAlgorithm object used to copy from
  */
  FEAlgorithm& operator=( const FEAlgorithm& p_feAlgorithm ) = delete;

  /*
   * Setup function that parses all the parameters needed by this class from the input file
   @param datafile - Reference to the input file
   @param sectionPath - a string representing the path to the section in the input file where the mesh parameters are specified
   @param sectionName - name of the section in the input file that holds all the information about the mesh parameters
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) override;
  
  /*
   * A function that sets up the problem parameters
   @param com - The communicator used to communicate between processors
   @param mesh - A pointer to the mesh on which the problem is defined
   @param fields - A pointer to the fields object that holds the fields involved in the problem
   @param prescribedConditions - A pointer to the object holding the prescribed conditions that define the problem
   @param nonlinearSolvers - A pointer to the object owning all the defined nonlinear solvers
   @param functions - A pointer to the object owning all the defined functions for this problem
   @param equations - A pointer to the object holding information about the weak form of the partial differential equation
   @param linearSystems - A pointer to the object holding information about the linear systems
   @param linearSolvers - A pointer to the object holding information about the linear solvers
  */
  void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com,
		   const NSMesh::Mesh* const p_mesh, const NSField::Fields* const p_fields,
		   const NSPrescribedCondition::PrescribedConditions* const p_prescribedConditions,
		   const NSNonlinearSolver::NonlinearSolvers* const p_nonlinearSolvers,
		   const NSInputFunction::Functions* const p_functions,
		   const NSEquation::Equations* const p_equations,
		   const NSLinearSystem::LinearSystems* const p_linearSystems,
		   const NSLinearSystem::LinearSolvers* const p_linearSolvers ) override;


  /*
   * Function that executes the finite element algorithm
   * The execution consists of the algorithm consists of
   * a) Updating the fields of the finite element problem
   * b) Assembling the linear system
   * c) Applying the strongly enforced conditions
   * d) Solving the linear system
   */
  void Execute() override;
  
private:
  /*
   * Function that sets information about all the fields of the problem
   * Fields of the problem include the unknowns, the test functions, the coordinates, and the input fields
   @param fields - A pointer to the fields object that holds the fields involved in the problem
  */
  void _initializeProblemFields( const NSField::Fields* const p_fields );

  /*
   * Function that sets the weak form of the problem
   @param equations - A pointer to the object holding information about the weak form of the partial differential equation
   @param ndim - Dimension of the problem
  */
  void _initializeWeakForm( const NSEquation::Equations* const p_equations, int p_ndim );

  /*
   * Function that sets information related to the linear system and linear solvers to be used to solve the problem
   @param com - The communicator used to communicate between processors
   @param linearSystems - A pointer to the object holding information about the linear systems
   @param linearSolvers - A pointer to the object holding information about the linear solvers
   @param fields - A pointer to the fields object that holds the fields involved in the problem
  */
  void _initializeLinearSystem( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com,
				const NSLinearSystem::LinearSystems* const p_linearSystems, const NSLinearSystem::LinearSolvers* const p_linearSolvers,
				const NSField::Fields* const p_fields ) override;

  /*
   * Helper function in the creation and setup of the linear system
   @param linearSystems - A pointer to the object holding information about the linear systems
   @param linearSolvers - A pointer to the object holding information about the linear solvers
  */
  void _createLinearSystem( const NSLinearSystem::LinearSystems* const p_linearSystems, const NSLinearSystem::LinearSolvers* const p_linearSolvers );
  
private:
  std::vector<std::string> _unknownsName = {}; /**<Names of the unknowns in the weak form*/
  std::vector<std::string> _testFunctionsName = {}; /**<Names of the test functions in the weak form*/
  std::vector<std::string> _inputsName = {}; /**<Names of the input to the weak form (body loads, material parameters, etc)*/
  std::unique_ptr<NSEquation::NSWeakForm::WeakForms> _weakForms = {}; /**<Pointer to the object representing the list of weak form of the pde problem*/
  
  std::map<std::string, NSField::Field<double>* > _nameToFieldMap = {}; /**<Map from the name of the fields involved in the finite element problem to a pointer to the field*/
  std::map<std::string, int> _nameToIdMap = {}; /**<Map from the name of a field involved in the fe problem to its id with respect to the total number of fields of the problem*/
  std::map<std::string, NSTypes::FieldOrder> _testFunctionNameToOrderMap = {}; /**<Map from the name of the test function to its order*/
  std::map<std::string, std::unique_ptr<NSField::Field<int> > > _nameToFieldDofIdsMap = {}; /**<Map from the name of the fields involved in the finite element problem to a pointer to the field holding the degrees of freedom ids corresponding to each field*/
  std::map<std::string, std::unique_ptr<NSField::FunctionSpace> > _nameToFieldDofsIdFs = {}; /**<Map from the name of a function space to a pointer to a function space object associated with the fields holding the ids of the dofs of the unknowns of the problem*/
  std::unique_ptr<NSFE::FEAssemblyOperator> _assemblyOperator = {}; /**<A pointer to the object responsible for assembling the lhs/rhs*/
};

/*
* Create function used by the factory to create this class
@output A pde algorithm whose static type is PDEAlgorithm and whose dynamic type is FEAlgorithm
*/
PDEAlgorithm* createFEAlgorithm();

namespace {
  static bool register_fe_algorithm = PDEAlgorithm::factory_pde_algorithm_Type::Register( "FE", &createFEAlgorithm );
}
    
} // namespace NSProblem
} // namespace NSSimulation
    
#endif
