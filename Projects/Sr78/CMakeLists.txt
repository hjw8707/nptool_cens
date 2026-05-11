cmake_minimum_required (VERSION 3.10)

project(Sr78)
# Setting the policy to match Cmake version
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
# include the default NPAnalysis cmake file
include("../../NPLib/ressources/CMake/NPAnalysis.cmake")

target_link_libraries (NPAnalysis -L${NPLIB}/lib -lNPCACAO -lNPTina -lNPDali2 -lNPGRAPE)
