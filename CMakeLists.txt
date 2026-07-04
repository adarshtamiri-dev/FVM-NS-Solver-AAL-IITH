cmake_minimum_required(VERSION 3.10)

project(Meshing_Phase)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_package(CGNS REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)


set(SOURCES 

main.cpp

preproc/PreProcessor.cpp
preproc/CGNSReader.cpp
preproc/GenerateFaceData.cpp
preproc/ComputeProperties.cpp
preproc/VerifyMesh.cpp

)

add_executable(LAB_SOLVERv1 ${SOURCES})


target_include_directories(LAB_SOLVERv1 PRIVATE /opt/opt/homebrew/include ${CGNS_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} .)
target_link_directories(LAB_SOLVERv1 PRIVATE /opt/homebrew/lib)

target_link_libraries(LAB_SOLVERv1 cgns hdf5 z ${CGNS_LIBRARIES} ${HDF5_LIBRARIES}

)