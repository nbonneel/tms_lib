cmake_minimum_required(VERSION 3.18)

project(tms_lib_project LANGUAGES CXX)

# Options
option(TMS_USE_CUDA "Enable CUDA code and define TMS_USE_CUDA" OFF)
option(TMS_USE_OPENMP "Enable OpenMP when available" ON)

# CUDA architecture list. For RTX 2080: 75.
# Override with: -DTMS_CUDA_ARCHITECTURES=75 or e.g. "75;86;89"
set(TMS_CUDA_ARCHITECTURES "75" CACHE STRING "CUDA architectures to compile for")

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Core library: CPU code
add_library(tms_lib STATIC
    tms_lib.cpp
    tms_lib.h
    SobolHardcodedTables.h
)

target_include_directories(tms_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# OpenMP, enabled by default but optional
if(TMS_USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(tms_lib PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(tms_lib PUBLIC TMS_USE_OPENMP)
        message(STATUS "OpenMP enabled")
    else()
        message(WARNING "TMS_USE_OPENMP=ON, but OpenMP was not found. Building without OpenMP.")
    endif()
else()
    message(STATUS "OpenMP disabled")
endif()

# Optional CUDA code
if(TMS_USE_CUDA)
    enable_language(CUDA)

    target_sources(tms_lib PRIVATE
        tms_lib_cuda_kernels.cu
    )

    target_compile_definitions(tms_lib PUBLIC TMS_USE_CUDA)

    set_target_properties(tms_lib PROPERTIES
        CUDA_STANDARD 14
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF
        CUDA_ARCHITECTURES "${TMS_CUDA_ARCHITECTURES}"
    )

    message(STATUS "CUDA enabled")
    message(STATUS "CUDA architectures: ${TMS_CUDA_ARCHITECTURES}")
else()
    message(STATUS "CUDA disabled")
endif()

# Demo / executable using main.cpp
add_executable(tms_demo
    main.cpp
)

target_link_libraries(tms_demo PRIVATE
    tms_lib
)

# Warnings / optimization niceties
if(MSVC)
    target_compile_options(tms_lib PRIVATE /W3)
    target_compile_options(tms_demo PRIVATE /W3)
else()
    target_compile_options(tms_lib PRIVATE -Wall -Wextra)
    target_compile_options(tms_demo PRIVATE -Wall -Wextra)
endif()
