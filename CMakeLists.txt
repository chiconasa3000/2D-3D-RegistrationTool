# CPU Implementation
set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 2.6)

include(${CMAKE_SOURCE_DIR}/IJMacros.txt)

project(CPUitkMultiImageRegistration)

#The following lines are required to use Dart
#enable_testing()
#include(Dart)

set_source_files_properties(*.txx PROPERTIES LANGUAGE CXX)

#Archivos ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#set(Required_Packages ITK)

#FOREACH(Package ${Required_Packages})
#  LOADPACKAGE(${Package})
##ENDFOREACH(Package)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
#set(Libraries ITKAlgorithms ITKBasicFilters ITKCommon ITKIO)
set(Libraries ${ITK_LIBRARIES})

file(GLOB TXX_FILES *.txx)

#set(TXX_DIR "/home/cris/Gpu-Midas-Journal-800")
#set(FILES_TXX ${TXX_DIR}/itkNormalizedGradientCorrelationImageToImageMetric.txx)

set(GPU_DIR "/Users/student/BasePrograms/InsightToolkit-4.13.0/Modules/Filtering/GPUImageFilterBase/include")
set(FILES_GPU ${GPU_DIR}/itkGPUNeighborhoodOperatorImageFilter.hxx ${GPU_DIR}/itkGPUNeighborhoodOperatorImageFilter.h)

set(HeaderFiles
itkGradientDifferenceMultiImageToImageMetric.h
itkGradientDifferenceSingleImageToImageMetric.h
itkGradientDifferenceSingleImageToImageMetric.txx
itkMeanSquaresMultiImageToImageMetric.h
itkMultiImageRegistrationMacro.h
itkMultiImageToImageMetric.h
itkMultiImageToImageMetric.txx
itkMultiImageToImageRegistrationMethod.h
itkMultiImageToImageRegistrationMethod.txx
itkMultiResolutionMultiImageToImageRegistrationMethod.h
itkMultiResolutionMultiImageToImageRegistrationMethod.txx
itkNormalizedGradientCorrelationImageToImageMetric.h
itkNormalizedGradientCorrelationImageToImageMetric.txx
itkNormalizedGradientCorrelationMultiImageToImageMetric.h
itkPatchedRayCastInterpolateImageFunction.h
itkPatchedRayCastInterpolateImageFunction.txx
itkPatternIntensityImageToImageMetric.h
itkPatternIntensityImageToImageMetric.txx
itkPatternIntensityMultiImageToImageMetric.h
itkPatternIntensityMultiImageToImageMetric.txx
)

#Add project's executables
#add_executable(ImageCompare ImageCompare.cxx)
#target_link_libraries(ImageCompare ${Libraries})

#add_executable(MultiImageSearch Examples/MultiImageSearch.cxx)
#target_link_libraries(MultiImageSearch ${Libraries})

add_executable(ReadCmd Examples/readCmd.cpp)
target_link_libraries(ReadCmd ${Libraries})

add_executable(MultiImageRegistration Examples/MultiImageRegistration.cxx  ${FILES_GPU} ${HeaderFiles})
target_link_libraries(MultiImageRegistration ${Libraries})


#install(FILES ${HeaderFiles} DESTINATION include COMPONENT Development)

#install(TARGETS MultiImageSearch DESTINATION bin COMPONENT Runtime)
install(TARGETS MultiImageRegistration DESTINATION bin COMPONENT Runtime)
#install(TARGETS ReadCmd DESTINATION bin COMPONENT Runtime)
