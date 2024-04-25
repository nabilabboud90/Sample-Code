#ifndef SOLUTIONTOFIELDDATAHELPER_HPP_
#define SOLUTIONTOFIELDDATAHELPER_HPP_

#include <string>
#include <map>
#include <memory>

namespace NSMesh { class Mesh; }
namespace NSField { template<typename T> class Field; }
namespace NSLinearSystem { class LinearSystem; }

namespace NSSimulation { namespace NSProblem { namespace NSFE {
      
void SetSolutionVectorFromFieldsData( const NSMesh::Mesh* p_mesh, const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
				      const std::map<std::string, std::unique_ptr<NSField::Field<int> > >& p_nameToFieldDofIdsMap,
				      NSLinearSystem::LinearSystem* const p_linearSystem );
      
void UpdateFieldDataFromSolutionIncrement( const NSMesh::Mesh* const p_mesh, const std::map<std::string, NSField::Field<double>* >& p_nameToFieldMap,
					   const std::map<std::string, std::unique_ptr<NSField::Field<int> > >& p_nameToFieldDofIdsMap,
					   NSLinearSystem::LinearSystem* const p_linearSystem );
      
} // namespace NSFE
} // namespace NSProblem
} // namespace NSSimulation

#endif
