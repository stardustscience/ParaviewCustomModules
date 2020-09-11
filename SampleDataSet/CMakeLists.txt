cmake_minimum_required(VERSION 3.9)
project(SampleDataset)

find_package(VTK)
include(${VTK_USE_FILE})

add_library(SampleDataset SHARED sampledataset.cpp point.cpp)
set_target_properties(SampleDataset PROPERTIES CXX_STANDARD 11 SUFFIX ".so")
install(TARGETS SampleDataset DESTINATION ParaviewPythonModules)
install(FILES SampleDataSet.py DESTINATION ParaviewPythonModules)

