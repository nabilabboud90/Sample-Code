#ifndef PDEPROBLEM_HPP_
#define PDEPROBLEM_HPP_

#include "Simulation/Problem/Problem.hpp"
#include "Mesh/Mesh.hpp"
#include "Mesh/STKMesh/STKMesh.hpp"
#include "Mesh/MeshQuery/MeshQuery.hpp"
#include "Mesh/Entity/Bucket.hpp"
#include "Field/Fields.hpp"
#include "PrescribedCondition/PrescribedConditions.hpp"
#include "NonlinearSolver/NonlinearSolvers.hpp"
#include "Output/Outputs.hpp"
#include "PDEAlgorithms/PDEAlgorithm.hpp"
#include "InputFunction/Functions.hpp"
#include "LinearSystem/LinearSolvers.hpp"
#include "LinearSystem/LinearSystems.hpp"
#include "Equation/Equations.hpp"

namespace NSSimulation { namespace NSProblem {

class PDEProblem : public Problem {
public:
  /*
   * Default constructor
   */
  explicit PDEProblem();

  /*
   * Copy constructor not implemented
   @param pdeProblem - PDEProblem object used to copy from
   */
  PDEProblem( const PDEProblem& p_pdeProblem ) = delete;

  /*
   * Equal operator not implemented
   @param pdeProblem - PDEProblem object used to copy from
  */
  PDEProblem& operator=( const PDEProblem& p_pdeProblem ) = delete;

  /*
   * Setup function that parses all the parameters needed by this class from the input file
   @param datafile - reference to the input file
   @param sectionPath - a string representing the path to the section in the input file where the simulation parameters are specified
   @param sectionName - name of the section in the input file that holds all the information about the simulation parameters
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) override;
  
  /*
   * Initialization function
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) override;

  /*
   * Function that executes the pde problem
   */
  void Execute() override;
  
private:
  std::unique_ptr<NSMesh::Mesh> _mesh = {}; /**<Pointer to the mesh object*/
  std::unique_ptr<NSField::Fields> _fields = {}; /**<Pointer to the fields object*/
  std::unique_ptr<NSInputFunction::Functions> _functions = {}; /**<A pointer to an object holding function expressions specified in the input file*/
  std::unique_ptr<NSPrescribedCondition::PrescribedConditions> _prescribedConditions = {}; /**<Pointer to the prescribed conditions (ic, bc, etc) object*/
  std::unique_ptr<NSNonlinearSolver::NonlinearSolvers> _nonlinearSolvers = {}; /**<Pointer to the nonlinear solvers objects*/
  std::unique_ptr<NSOutput::Outputs> _outputs = {}; /**<Pointer to the outputs object*/
  std::unique_ptr<NSLinearSystem::LinearSolvers> _linearSolvers = {}; /**<Pointer to the object holding data about the linear solvers defined for this problem*/
  std::unique_ptr<NSLinearSystem::LinearSystems> _linearSystems = {}; /**<Pointer to the object holding data about the linear systems defined for this problem*/
  std::unique_ptr<NSEquation::Equations> _equations = {}; /**<Pointer to the object holding information about this equations of this problem*/
  
  std::vector<std::string> _algorithmNames = {}; /**<A vector representing the names of the algoithms*/
  std::map<std::string, std::unique_ptr<PDEAlgorithm> > _nameToAlgorithmMap = {}; /**<Map from the name of the algorithm to the pointer to the object representing the algorithm*/
};

/*
* Create function used by the factory to create this class
@output A problem whose static type is Problem and whose dynamic type is PDEProblem
*/
Problem* createPDEProblem();

namespace {
  static bool register_PDEProblem = Problem::factory_problem_Type::Register( "PDEProblem", &createPDEProblem );
}

      
} // namespace NSProblem
} // namespace NSSimulation

#endif
