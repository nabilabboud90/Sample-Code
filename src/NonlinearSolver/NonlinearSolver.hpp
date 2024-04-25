#ifndef NONLINEARSOLVER_HPP_
#define NONLINEARSOLVER_HPP_

#include <string>
#include <memory>
#include <functional>

namespace NSFactory { template<typename Product, typename Identifier, typename... Argument> class Factory; }
namespace NSParser { namespace NSInput { class Datafile; } }
namespace NSCommunicator { namespace NSMPIComWrapper { class MPIComWrapper; } }

namespace NSNonlinearSolver {

class NonlinearSolver {
public:
  typedef NSFactory::Factory<NonlinearSolver*, std::string, void> factory_nonlinearsolver_Type; /* *< Factory to chose among the different nonlinear solvers supported in the code*/
  
  /*
   * Default constructor
   */
  explicit NonlinearSolver();

  /*
   * Copy constructor not implemented
   @param nsd - NonlinearSolver class to copy from
  */
  NonlinearSolver( const NonlinearSolver& p_nsd ) = delete;

  /*
   * Equal operator not implemented
   @param nsd - NonlinearSolver class to set equal to
  */
  NonlinearSolver& operator=( const NonlinearSolver& p_nsd ) = delete;

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - Reference to the input file
   @param sectionPath - A string representing the path to the section in the input file
   @param sectionName - Name of the section in the input file
  */
  virtual void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Initialization function
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  virtual void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com );

  /*
   * Function to run one newton iteration
   @param linearize - Pointer to the function that will assemble and update the linearized system of equations representing the nonlinear equations to be solved
   @param finalize - Pointer to the function that will do any post assembly work needed by the linear system representing the nonlinear problem solved
   @param solve - Pointer to the function that will call a solve method on the linear system representing the nonlinear problem
   @param print - Pointer to the function that prints the update about the current newton iteration
   @param verbose - Boolean to allow whether to print or not the norms of the residual
   @output Boolean specifying whether the nonlinear problem has converged
  */
  virtual bool Iterate( std::function<void()> p_linearize, std::function<std::vector<double>()> p_finalize, std::function<void()> p_solve,
			std::function<void( double, double, int, double )> p_print, bool verbose ) const = 0;
  
protected:
  std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper> _com = {}; /**<Communicator responsible for doing any communication on the mesh*/
  std::string _name = {}; /**<String holding the name of the object which will be used to reference it in other places of the code*/
  int _maxIterations = {}; /**<An integer representing the maximum number of iterations allowed for convergence*/
  double _relativeTolerance = {}; /**<Double representing the relative tolerance used to check for convergence*/
  double _absoluteTolerance = {}; /**<Double representing the absolute tolerance used to check for convergence*/
};
  
} // namespace NSNonlinearSolver

#endif
