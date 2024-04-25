#include<iostream>
#include<mpi.h>

#include "src/Communicator/MPI/MPIComWrapper.hpp"
#include "src/Utility/Timer.hpp"
#include "Parser/InputParser/Datafile.hpp"
#include "Simulation/Simulation.hpp"

#include <Kokkos_Core.hpp>

int main( int p_argc, char* p_argv[] ) {
  // Check that the code only takes one argument
  if( p_argc != 2 ) {
    throw std::invalid_argument( "The code only takes one argument corresponding to the input file" );
  }

  // Initialize MPI
  MPI_Init( nullptr, nullptr );
  
  // Initialize Kokkos
  Kokkos::initialize( p_argc, p_argv );
  {    
    // Define the timer
    NSTimer::Timer timer;
    timer.Start( "Total Simulation Time" );
    
    // Create an mpi communicator wrapper
    std::shared_ptr<NSCommunicator::NSMPIComWrapper::MPIComWrapper> mpi_com( new NSCommunicator::NSMPIComWrapper::MPIComWrapper );
    mpi_com->SetMPICom( MPI_COMM_WORLD );
    
    // Parse the datafile
    std::string datafile_name( p_argv[1] );
    std::unique_ptr<NSParser::NSInput::Datafile> datafile = std::make_unique<NSParser::NSInput::Datafile>();
    datafile->Read( datafile_name );
    
    // Create the problem
    NSSimulation::Simulation simulation;

    // Setup the problem
    simulation.Setup( *datafile, "", "SIMULATION" );

    // Initialize the probelm
    simulation.Initialize( mpi_com );

    // Execute the simulation
    simulation.Execute();
    
    timer.End( "Total Simulation Time" );
  }
  
  // Finalize Kokkos
  Kokkos::finalize();
  
  // Finalize MPI
  MPI_Finalize();

  return 1;
}
