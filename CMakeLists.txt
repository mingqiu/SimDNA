cmake_minimum_required(VERSION 3.3)
project(SimDNA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp lib/Vector3D.h include/DataStructure.h src/DataStructure.cpp include/Origami.h src/Origami.cpp include/Configure.h src/Test.cpp src/Origami_functions.cpp src/ToOpenMM.cpp src/Utility.cpp include/Utility.h src/TrashCan.cpp)
add_executable(SimDNA ${SOURCE_FILES})