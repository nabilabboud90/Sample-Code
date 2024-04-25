#ifndef NEWTONSOLVER_HPP_
#define NEWTONSOLVER_HPP_

#include <memory>
#include <functional>

#include "NonlinearSolver/NonlinearSolver.hpp"
#include "Factory/Factory.hpp"

namespace NSParser { namespace NSInput { class Datafile; } }
namespace NSCommunicator { namespace NSMPIComWrapper { class MPIComWrapper; } }

namespace NSNonlinearSolver {

class NewtonSolver : public NonlinearSolver {
public:
  /*
   * Default constructor
   */
  explicit NewtonSolver();

  /*
   * Copy constructor not implemented
   @param ns - NewtonSolver class to copy from
  */
  NewtonSolver( const NewtonSolver& p_ns ) = delete;

  /*
   * Equal operator not implemented
   @param ns - LienarSolver class to set equal to
  */
  NewtonSolver& operator=( const NewtonSolver& p_ns ) = delete;

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - Reference to the input file
   @param sectionPath - A string representing the path to the section in the input file where the newton solver parameters are specified
   @param sectionName - Name of the section in the input file that holds all the information about the newton solver parameters
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Initialization function
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) override;

  /*
   * Function to run one newton iteration
   @param linearize - Pointer to the function that will assemble and update the linearized system of equations representing the nonlinear equations to be solved
   @param finalize - Pointer to the function that will do any post assembly work needed by the linear system representing the nonlinear problem solved
   @param solve - Pointer to the function that will call a solve method on the linear system representing the nonlinear problem
   @param print - Pointer to the function that prints the update about the current newton iteration
   @param verbose - Boolean to allow whether to print or not the norms of the residual
   @output Boolean specifying whether the nonlinear problem has converged
  */
  bool Iterate( std::function<void()> p_linearize, std::function<std::vector<double>()> p_finalize, std::function<void()> p_solve,
		std::function<void( double, double, int, double )> p_print, bool verbose ) const override;

private:
  double _etaMaximum = {}; /**<Tolerance used to update the linear solver tolerances based on the residual norm*/
  double _gamma = {}; /**<Tolerance used to update the linear solver tolerances based on the residual norm*/
};

/*
* Create function used by the factory to create this class
@output A nonlinear solver whose static type is NonlinearSolver and whose dynamic type is NewtonSolver
*/
NonlinearSolver* createNewtonSolver();

namespace {
  static bool register_newton_solver = NonlinearSolver::factory_nonlinearsolver_Type::Register( "NEWTON", &createNewtonSolver );
}
  
}

#endif
