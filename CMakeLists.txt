cmake_minimum_required (VERSION 2.6)

set(CMAKE_BUILD_TYPE "Debug")

project (DaFFEM)
enable_language (Fortran)

if ( ${USE_MUMPS} )
   message(" ")
   message(STATUS "Compiling with MUMPS")
   message(STATUS "Home directory: " ${MUMPS_DIR})
   add_definitions(-DMUMPS)
   message(STATUS "include directory: " ${MUMPS_DIR}/include) 
   include_directories( "${MUMPS_DIR}/include")
   message(STATUS "link directories: " ${MUMPS_DIR}/lib "; " ${MUMPS_DIR}/libseq)
   link_directories( "${MUMPS_DIR}/lib" "${MUMPS_DIR}/libseq")
   message(" ")
 endif()
 link_directories( "/usr/local/lib")
 add_definitions("-fbounds-check")
 add_definitions("-cpp")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${add_flag}")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}" "./mod")
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# file(GLOB SOURCES "src/*.f90" "src/SparseMatrix/*.f90" "src/Embedded/*.f90" "src/Schemes/*.f90" "src/SystemResolution/*.f90" "src/SystemResolution/*.F90" "src/GradientReconstruction/*.f90")
file(GLOB SOURCES "src/*.f90" "src/*/*.f90" "src/*/*/*.f90" "src/*/*/*/*.f90")

add_executable(DaFFEM ${SOURCES})
if ( USE_MUMPS )
  target_link_libraries(DaFFEM dmumps mumps_common pord mpiseq blas pthread)
endif()
