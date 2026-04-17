cmake_minimum_required(VERSION 3.16)

# Set macOS deployment target before project() to fix gfortran linker issues
if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0" CACHE STRING "macOS deployment target")
    set(CMAKE_OSX_SYSROOT "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk")
endif()

project(loglikelihood LANGUAGES Fortran)

find_package(MPI REQUIRED)
if(MPI_Fortran_FOUND)
    enable_language(Fortran)
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
endif()

# Main loglikelihood executable
add_executable(loglikelihood
    compute_loglikelihood.f90
    main_loglikelihood.f90
)

# Batch loglikelihood executable
add_executable(batch_loglikelihood
    compute_loglikelihood.f90
    batch_loglikelihood.f90
)

target_link_libraries(loglikelihood ${MPI_Fortran_LIBRARIES})
target_link_libraries(batch_loglikelihood ${MPI_Fortran_LIBRARIES})
