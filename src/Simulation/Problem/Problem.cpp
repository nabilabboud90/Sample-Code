#include "Problem.hpp"
#include "Parser/InputParser/Datafile.hpp"

namespace NSSimulation {

Problem::Problem() {
}


void
Problem::Setup( [[maybe_unused]]const NSParser::NSInput::Datafile& p_datafile, [[maybe_unused]]const std::string& p_sectionPath, [[maybe_unused]]const std::string& p_sectionName ) {
}


void
Problem::Initialize( const std::shared_ptr<const NSCommunicator::NSMPIComWrapper::MPIComWrapper>& p_com ) {
  _com = p_com;
}

} // namespace NSSimulation
