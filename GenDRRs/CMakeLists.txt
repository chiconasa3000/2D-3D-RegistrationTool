cmake_minimum_required(VERSION 2.8)

project(GenDRRs)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#generacion del ejeutable
add_executable(genVirtualImage genVirtualImage.cxx)

#Vinvulo del ejecutable con la libreria ITK
target_link_libraries(genVirtualImage ${ITK_LIBRARIES})
