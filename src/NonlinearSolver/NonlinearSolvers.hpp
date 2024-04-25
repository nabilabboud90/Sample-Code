#ifndef NONLINEARSOLVERS_HPP_
#define NONLINEARSOLVERS_HPP_

#include<memory>
#include<vector>
#include<map>

#include "NonlinearSolver/NonlinearSolver.hpp"

namespace NSParser { namespace NSInput { class Datafile; } }
namespace NSCommunicator { namespace NSMPIComWrapper { class MPIComWrapper; } }

namespace NSNonlinearSolver {

class NonlinearSolvers {
public:
  /*
   * Default constructor
   */
  explicit NonlinearSolvers();

  /*
   * Copy constructor not implemented
   @param nsd - NonlinearSolvers class to copy from
  */
  NonlinearSolvers( const NonlinearSolvers& p_nsd ) = delete;

  /*
   * Equal operator not implemented
   @param nsd - NonlinearSolvers class to set equal to
  */
  NonlinearSolvers& operator=( const NonlinearSolvers& p_nsd ) = delete;

  /*
   * Setup function that mainly parses all the parameters needed by this class from the input file
   @param datafile - Reference to the input file
   @param sectionPath - A string representing the path to the section in the input file
   @param sectionName - Name of the section in the input file
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Initialization function
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com );

  /*
   * Get a newton solver of a specified name
   @param name - Name of the newton solver to be returned
   @output Pointer to the newton solver of a specified name
  */
  const NSNonlinearSolver::NonlinearSolver* GetNonlinearSolver( const std::string& p_name ) const;

private:
  std::vector<std::string> _nonlinearSolversName = {}; /**<Vector holding the names of the nonlinear solvers*/
  std::map<std::string, std::unique_ptr<NSNonlinearSolver::NonlinearSolver> > _nameToNonlinearSolverMap; /**<Map from the nonlinear solver name to a pointer to it*/
};

} // namespace NonlinearSolver

#endif
