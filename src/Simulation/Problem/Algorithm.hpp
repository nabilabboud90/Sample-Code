#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_

#include <string>

namespace NSSimulation { namespace NSProblem {

class Algorithm {
public:
  /*
   * Default constructor
   */
  explicit Algorithm();

  /*
   * Copy constructor not implemented
   @param algorithm - Algorithm object used to copy from
   */
  Algorithm( const Algorithm& p_algorithm ) = delete;

  /*
   * Equal operator not implemented
   @param algorithm - Algorithm object used to copy from
  */
  Algorithm& operator=( const Algorithm& p_algorithm ) = delete;

private:
  std::string _name = {}; /**<Name of the algorithm used to reference it in the code*/
};

} // namespace NSProblem
} // namespace NSSimulation

#endif 
