#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_

#include <string>
#include <memory>

namespace NSCommunicator { namespace NSMPIComWrapper { class MPIComWrapper; } }
namespace NSParser { namespace NSInput { class Datafile; } }
namespace NSFactory { template<typename Product, typename Identifier, typename... Argument> class Factory; }

namespace NSSimulation {

class Problem {
public:
  typedef NSFactory::Factory<Problem*, std::string, void> factory_problem_Type; /* *< Factory to chose among the different problem types supported in the code */
  
  /*
   * Default constructor
   */
  explicit Problem();

  /*
   * Copy constructor not implemented
   @param problem - Problem object used to copy from
   */
  Problem( const Problem& p_problem ) = delete;

  /*
   * Equal operator not implemented
   @param problem - Problem object used to copy from
  */
  Problem& operator=( const Problem& p_problem ) = delete;

  /*
   * Setup function that parses all the parameters needed by this class from the input file
   @param datafile - reference to the input file
   @param sectionPath - a string representing the path to the section in the input file where the simulation parameters are specified
   @param sectionName - name of the section in the input file that holds all the information about the simulation parameters
  */
  virtual void Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName );
  
  /*
   * Initialization function
   @param com - The communicator that will be invoked whenever any parallel communication need to be done
  */
  virtual void Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com );

  /*
   * Function that executes the problem
   */
  virtual void Execute() = 0;
  
protected:
  std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper> _com = {}; /**<Communicator responsible for doing any communication between processors*/
};

} // namespace NSSimulation

#endif
