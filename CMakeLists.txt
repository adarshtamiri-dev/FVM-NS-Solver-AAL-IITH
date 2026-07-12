cmake_minimum_required(VERSION 3.10)

project(Meshing_Phase)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CGNS_ROOT ${CMAKE_SOURCE_DIR}/preproc/cgns)
set(CMAKE_PREFIX_PATH ${CGNS_PREFIX_PATH} ${CGNS_ROOT})

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
add_executable(LAB_SOLVERv1 main.cpp)
target_include_directories(LAB_SOLVERv1 PRIVATE ${CGNS_INCLUDE_DIRS})
target_link_libraries(LAB_SOLVERv1 PRIVATE ${CGNS_LIBRARIES})
