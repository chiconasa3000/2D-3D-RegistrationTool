cmake_minimum_required(VERSION 2.8)

project(GenDRRs)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#recognize txx and cxx files
set_source_files_properties(*.txx PROPERTIES LANGUAGE CXX)

#Libraries ITK
set(Libraries ${ITK_LIBRARIES})

#consider TXX files
file(GLOB TXX_FILES *.txx)

#setting the headers
set(HeaderFiles 
../itkPatchedRayCastInterpolateImageFunction.h
../itkPatchedRayCastInterpolateImageFunction.txx)

#generacion del ejeutable
add_executable(genVirtualImage genVirtualImage.cxx)

#Vinvulo del ejecutable con la libreria ITK
target_link_libraries(genVirtualImage ${ITK_LIBRARIES})





