#ifndef FEPRESCRIBEDCONDITIONSHELPER_HPP_
#define FEPRESCRIBEDCONDITIONSHELPER_HPP_

#include <string>
#include <map>
#include <memory>

namespace NSMesh { class Mesh; }
namespace NSField { template<typename T> class Field; }
namespace NSPrescribedCondition { class PrescribedCondition; }
namespace NSLinearSystem { class LinearSystem; }

namespace NSSimulation { namespace NSProblem { namespace NSFE {

void ApplyBoundaryConditionsOnFields( const NSMesh::Mesh* p_mesh, const NSField::Field<double>* p_coordinates,
				      const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
				      const std::map<std::string, const NSPrescribedCondition::PrescribedCondition*>& p_nameToBoundaryConditionMap,
				      double p_time );

void ApplyInitialConditionsOnFields( const NSMesh::Mesh* p_mesh, const NSField::Field<double>* p_coordinates,
				     const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
				     const std::map<std::string, const NSPrescribedCondition::PrescribedCondition*>& p_nameToInitialConditionMap,
				     double p_time );

void ApplyBoundaryConditionsOnLinearSystem( const NSMesh::Mesh* p_mesh, const std::map<std::string, std::unique_ptr<NSField::Field<int> > >& p_nameToFieldDofIdsMap,
					    const std::map<std::string, const NSPrescribedCondition::PrescribedCondition*>& p_nameToBoundaryConditionMap,
					    NSLinearSystem::LinearSystem* const p_linearSystem );
      
} // namespace NSFE
} // namespace NSProblem
} // namespace NSSimulation

#endif
