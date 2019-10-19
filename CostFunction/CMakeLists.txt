cmake_minimum_required (VERSION 3.12.0)

project(MeanSquaresMetric)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


set(HeaderFiles 
	itkNormalizedGradientCorrelationImageToImageMetric.h
	itkNormalizedGradientCorrelationImageToImageMetric.txx
	utils.h
	itkPatchedRayCastInterpolateImageFunction.h	
)

add_executable(MeanSquaresMetric MeanSquaresImageMetric1.cxx utils.cpp)
target_link_libraries(MeanSquaresMetric ${Glue} ${VTK_LIBRARIES} ${ITK_LIBRARIES})
