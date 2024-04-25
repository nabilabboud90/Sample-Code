#include "Simulation.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "Problem/PDEProblem/PDEProblem.hpp"
#include "Factory/Factory.hpp"

namespace NSSimulation {

Simulation::Simulation() {
}


void
Simulation::Setup( const NSParser::NSInput::Datafile& p_datafile, const std::string& p_sectionPath, const std::string& p_sectionName ) {
  // Read the list of problems
  p_datafile.Get( p_sectionPath + "/" + p_sectionName + "/LIST", _problemsName );
  
  // Loop over each problem
  for( size_t p=0; p<_problemsName.size(); ++p ) {
    // Read the problem type
    std::string problem_type = p_datafile.Get<std::string>( p_sectionPath + "/" + p_sectionName + "/" + _problemsName[p] + "/TYPE" );
    
    // Read the datafile name that holds the information for this problem
    std::string datafile_name = p_datafile.Get<std::string>( p_sectionPath + "/" + p_sectionName + "/" + _problemsName[p] + "/DATAFILE" );
    
    // Create the datafile
    auto datafile = std::make_unique<NSParser::NSInput::Datafile>();
    
    // Read the datafile
    datafile->Read( datafile_name );
    
    // Create the problem
    std::unique_ptr<NSSimulation::Problem> problem;
    problem.reset( NSSimulation::Problem::factory_problem_Type::Create( problem_type ) );
    
    // Setup the problem
    problem->Setup( *datafile, "/" + _problemsName[p], "" );
    
    // Store the problem
    _nameToProblemMap.insert( std::make_pair( _problemsName[p], std::move( problem ) ) );
  }
}


void
Simulation::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  // Loop over each problem and initialize it
  for( auto& problem : _problemsName ) {
    // Initialize the problem
    _nameToProblemMap.find( problem )->second->Initialize( p_com );
  }
}


void
Simulation::Execute() {
  for( auto& problem : _nameToProblemMap ) {
    problem.second->Execute();
  }
}
  
} // namespace NSSimulation
