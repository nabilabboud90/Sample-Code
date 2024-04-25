# Specify the cmake version being used
cmake_minimum_required(VERSION 3.15.5)

# Name your project here
project(Code)

# Determine the operating system on which the code is being compiled
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
    set(MACOSX FALSE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(LINUX FALSE)
    set(MACOSX TRUE)
endif()
message (STATUS "Operating System: Linux = ${LINUX}")
message (STATUS "Operating System: Darwin= ${MACOSX}")

# Set the compilers for Code
MESSAGE("Setting and checking of compilers:")
SET(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${CMAKE_C_COMPILER} )

# Set the path to the trilinos include
include_directories(${Trilinos_Include_Dir})

# Set the path to the open mpi include
# Note i think because open mpi is installed as a shared library I only need to specify the path to the library in the system's environment variables
include_directories(${MPI_Include_Dir})

# Set the path to the boost include
# Note i think because boost is installed as a shared library I only need to specify the path to the library in the system's environment variables
include_directories(${Boost_Include_Dir})

# Silence some compiler generated warning about the boost library
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

# Set compiler flags
SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} )
string (REPLACE ";" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Change some settings if the build mode is debug
if(CMAKE_BUILD_TYPE MATCHES DEBUG)
message(STATUS "Removing flag ... ")
string (REPLACE "-O3" "-O0" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# Compile the source code
add_subdirectory(src)

# Add the executable
add_executable(Code main.cpp)

# Link the rest of the source code which was built as separate libraries to the main
target_link_libraries(Code combinedLib)

# Include the directories for the different subcomponents of the code
include_directories(${Code_Src_Dir}/src/)

# Link to the static libraries of trilinos
# I think that find_package will populate the Trilinos_LIBRARIES variable that i can use after to link to trilinos
find_package(Trilinos REQUIRED)
target_link_libraries(Code ${Trilinos_LIBRARIES})

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/Users/nabilabboud/Documents/TPL/lib/boost-1.68.0/include")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/Users/nabilabboud/Documents/TPL/lib/boost-1.68.0/lib")

# Link to the boost library
find_package(Boost REQUIRED)
target_link_libraries(Code ${Boost_LIBRARIES})

# Link to the lapack library
find_package(LAPACK REQUIRED)
target_link_libraries(Code ${LAPACK_LIBRARIES})

# Link to the doxygen library
#find_package(Doxygen REQUIRED)
#target_link_libraries(Code ${Doxygen_LIBRARIES})

# Link to the Vc library if it is available
if(${Enable_Vc})
message("\"Linking to Vc\"")
#find_package(Vc REQUIRED)
find_package(Vc ${Vc_FIND_VERSION} QUIET NO_MODULE PATHS $ENV{HOME} ${Vc_Lib_Dir})
target_link_libraries(Code ${Vc_LIBRARIES})
endif()
