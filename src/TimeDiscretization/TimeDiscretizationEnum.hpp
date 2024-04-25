#ifndef TIMEDISCRETIZATIONENUM_HPP_
#define TIMEDISCRETIZATIONENUM_HPP_

#include <string>

namespace NSTimeDiscretization {

enum Scheme { ts_QUASISTATIC, ts_BDF1, ts_BDF2, ts_FIRST = ts_QUASISTATIC, ts_LAST = ts_BDF2 };

NSTimeDiscretization::Scheme GetScheme( const std::string& p_scheme );
  
} // namespace NSTimediscretization

#endif
