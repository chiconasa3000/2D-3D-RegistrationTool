
#include <iostream>
#include "itkImage.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"

#include "itkSimilarity3DTransform.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkImageFileWriter.h"
#include "../itkPatchedRayCastInterpolateImageFunction.h"
#include "itkMatrix.h"
#include "gestorMatrixRot.h"
#include <string>
#include <itksys/SystemTools.hxx>
#include <fstream>
#include "utils.h"

const unsigned int Dimension = 3;
//constant for casting degrees into radians format of rotation projection
const double dtr = ( atan(1.0) * 4.0 ) / 180.0;

class GenVirtGen{
private:
	//Variables para la conversion al Versor de una rotacion Euler
	double ax,ay,az,angle;


	itk::TimeProbesCollectorBase timer;	//Time Record

	char *input_name = NULL;
	char *logFileName = NULL;
	
	std::ofstream logregistro;

	//tipo de pixel by default para las imagenes
	typedef short int InputPixelType;
	typedef unsigned char OutputPixelType;

	//typedef itk::Image<short int, Dimensions> FixedImageType;
	typedef itk::Image<InputPixelType, Dimension> MovingImageType;
	typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

	MovingImageType::Pointer image;

	typedef itk::ResampleImageFilter<MovingImageType, MovingImageType> FilterType;
	
	typedef itk::Similarity3DTransform< double > TransformType;
	TransformType::Pointer transform; 
		
	//Lectura de Parametros de transformacion para mostrarlo en la info
	TransformType::ParametersType similarityParameters;
	
	//Instance of the interpolator
	typedef itk::PatchedRayCastInterpolateImageFunction<MovingImageType, double> InterpolatorType;
	FilterType::Pointer filter;

	typedef InterpolatorType::InputPointType FocalPointType;
	
	InterpolatorType::Pointer interpolator;
	
	HelperRot helperRot;

	char *type_projection = NULL;	


public:
	GenVirtGen(char * typeProjection,char *logFileName);
	void setLogFile(char *logfile);	
	int readMovingImage(char *input_name);
	void initResampleFilter(int &dxx,int &dyy, float &im_sx, float &im_sy,bool customized_2DCX,float o2Dx, float o2Dy,float scd);
	void initTransform(float tx,float ty,float tz, float dcx, float dcy, float dcz,float rx, float ry, float rz,float sg);	
	void copyPropertiesInputImage(int &dxx, int &dyy, float &im_sx, float &im_sy);
	void initInterpolator(float focalPointx, float focalPointy, float focalPointz, float threshold);
	void writeResultImage(char *output_name);
	MovingImageType::Pointer returnResultImage(char *output_name);
};



