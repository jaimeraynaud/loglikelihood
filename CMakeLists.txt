
cmake_minimum_required(VERSION 3.16)

# Set macOS deployment target before project() to fix gfortran linker issues
if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0" CACHE STRING "macOS deployment target")
    set(CMAKE_OSX_SYSROOT "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk")
endif()

project(loglikelihood LANGUAGES Fortran)

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
