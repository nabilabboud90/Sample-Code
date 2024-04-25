#ifndef DISCRETIZEINTIME_HPP_
#define DISCRETIZEINTIME_HPP_

#include <vector>

#include "TimeDiscretization/TimeDiscretizationEnum.hpp"

namespace NSLinearSystem { class LinearSystem; }

namespace NSSimulation { namespace NSProblem {

void DiscretizeInTime( NSLinearSystem::LinearSystem* const p_linearSystem, const NSTimeDiscretization::Scheme& p_scheme, const std::vector<double>& p_dts );

} // namespace NSProblem
} // namespace NSSimulation

#endif
