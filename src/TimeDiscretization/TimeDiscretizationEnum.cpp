#include "TimeDiscretizationEnum.hpp"
#include "Utility/Exceptions.hpp"

#include <map>

namespace NSTimeDiscretization {

NSTimeDiscretization::Scheme GetScheme( const std::string& p_scheme ) {
  std::map<std::string, Scheme> result;
  result.insert( std::make_pair( "QUASISTATIC", ts_QUASISTATIC ) );
  result.insert( std::make_pair( "BDF1", ts_BDF1 ) );
  result.insert( std::make_pair( "BDF2", ts_BDF2 ) );

  NSExceptions::InvalidArgument( result.find( p_scheme ) == result.end(), "The time discretization scheme asked for is not defined in the TimeDiscretization Scheme enum" );

  return result.find( p_scheme )->second;
}
  
} // namespace NSTimediscretization
