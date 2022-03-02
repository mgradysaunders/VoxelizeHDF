# CMake version.
cmake_minimum_required(VERSION 3.14)

# Project.
project(
    VoxelizeHDF
    VERSION 1.0.0 
    DESCRIPTION "DIRSIG5 HDF Voxelizer"
    LANGUAGES C CXX
    )

set(CMAKE_C_COMPILER gcc-11)
set(CMAKE_CXX_COMPILER g++-11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_NATIVE "Build and tune for native hardware?" OFF)

# Set release.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Include GNUInstallDirs.
include(GNUInstallDirs)
message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "INSTALL_INCLUDEDIR: ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "INSTALL_LIBDIR: ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "INSTALL_BINDIR: ${CMAKE_INSTALL_BINDIR}")

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(Precept REQUIRED)
find_package(HDF5 1.10.0 COMPONENTS C CXX HL REQUIRED)

add_executable(
    VoxelizeHDF
    main.cpp
    )
set_target_properties(
    VoxelizeHDF
    PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED YES
    )
target_include_directories(
    VoxelizeHDF 
  PUBLIC
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(
    VoxelizeHDF 
  PUBLIC
    Precept::Precept
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
)
target_link_libraries(VoxelizeHDF PUBLIC -ldl -static)

install(
    TARGETS ${PROJECT_NAME} 
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
install(FILES ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} DESTINATION ${CMAKE_INSTALL_LIBDIR})
