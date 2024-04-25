#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Simulation/Problem/Problem.hpp"

namespace NSCommunicator { namespace NSMPIComWrapper { class MPIComWrapper; } }

namespace NSSimulation {

class Simulation {
public:
  /*
   * Default constructor
   */
  explicit Simulation();

  /*
   * Copy constructor not implemented
   @param simulation - Simulation object used to copy from
   */
  Simulation( const Simulation& p_simulation ) = delete;

  /*
   * Equal operator not implemented
   @param simulation - Simulation object used to copy from
  */
  Simulation& operator=( const Simulation& p_simulation ) = delete;

  /*
   * Setup function that parses all the parameters needed by this class from the input file
   @param datafile - reference to the input file
   @param sectionPath - a string representing the path to the section in the input file where the simulation parameters are specified
   @param sectionName - name of the section in the input file that holds all the information about the simulation parameters
  */
  void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );

  /*
   * Initialization function
   @param mpi_com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com );

  /*
   * Function that executes the simulation
   */
  void Execute();
  
private:
  std::vector<std::string> _problemsName = {}; /**<Vector holding the name of the problems that make up the simulation*/
  std::map<std::string, std::unique_ptr<NSSimulation::Problem> > _nameToProblemMap = {}; /**<Map from the name of a problem to the pointer to the problem. A problem can be a finite element problem, a finite difference, finite volume, modal problem, etc*/
};
  
} // namespace NSSimulation

#endif
