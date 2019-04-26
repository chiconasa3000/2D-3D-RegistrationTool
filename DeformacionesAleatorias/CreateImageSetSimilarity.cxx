/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: ResampleImageFilter.cxx,v $
Language:  C++
Date:      $Date: 2006/05/14 12:12:52 $
Version:   $Revision: 1.32 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"


#include <itkSimilarity3DTransform.h>
#include <itkEuler3DTransform.h>
#include "itkBSplineDeformableTransform.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkTransformFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkFlipImageFilter.h"
#include <string>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <itksys/SystemTools.hxx>

#include "utils.h"

using namespace std;

void menu(){
	std::cout << "Deformador de Imagenes Medicas "<< std::endl;
	std::cout << "inputImageFile -v -rnd -rnd_sem -folderName -numImages -rx -ry -rz -t -sg -centerSpec -originSpec"<< std::endl; 
}

int main( int argc, char * argv[] )
{

	//Versor 3D Rotation
	double ax,ay,az,angle;

	//Converter From Degrees to Radians
	const double dtr = (atan(1.0) * 4.0)/180.0;
	const double rtd = (180.0)/(atan(1.0)*4.0);
	//inputImageFile
	char *inputImageFile = NULL;
	//folderName
	char *folderName = NULL;
	//Log File Name in order to save the output stream
	char *logfilename = NULL; 	
	//numberOfImages
	int numberOfImages = 0;
	//Random Mode
	bool random_mode = false;
	//Semilla Random
	bool semillaRandom = false;
	//Rotacion
	float rx = 0.;
	float ry = 0.;
	float rz = 0.;	
	//Translation
	float tx = 0.;
	float ty = 0.;
	float tz = 0.;
	//Scale
	float sg = 1.;
	//Origin Inside
	float ox = 0.;
	float oy = 0.;
	float oz = 0.;
	
	//Origin by user
	bool originSpec = false;
	//Center in the corner
	bool centerSpec = false;
	
	bool ok; //sanity of parameters
	
	bool verbose = false;
	

	while(argc > 1){
		ok = false;

		if((ok ==false) && (strcmp(argv[1], "-v") == 0)){
			argc--; argv++;
			ok == true;
			verbose = true;
		}
		if((ok ==false) && (strcmp(argv[1], "-centerSpec") == 0)){
			argc--; argv++;
			ok == true;
			centerSpec = true;
		}
		
		if((ok ==false) && (strcmp(argv[1], "-rnd") == 0)){
			argc--; argv++;
			ok == true;
			random_mode = true;
		}

		if((ok ==false) && (strcmp(argv[1], "-rnd_sem") == 0)){
			argc--; argv++;
			ok == true;
			semillaRandom = true;
		}

		if((ok ==false) && (strcmp(argv[1], "-folderName") == 0)){
			argc--; argv++;
			ok == true;
			folderName = argv[1];
			argc--; argv++;
		}

		if((ok ==false) && (strcmp(argv[1], "-numImages") == 0)){
			argc--; argv++;
			ok == true;
			numberOfImages = atoi(argv[1]);
			argc--; argv++;
		}

		if((ok ==false) && (strcmp(argv[1], "-rx") == 0)){
			argc--; argv++;
			ok == true;
			rx = atof(argv[1]);
			argc--; argv++;
		}
		if((ok ==false) && (strcmp(argv[1], "-ry") == 0)){
			argc--; argv++;
			ok == true;
			ry = atof(argv[1]);
			argc--; argv++;
		}
		if((ok ==false) && (strcmp(argv[1], "-rz") == 0)){
			argc--; argv++;
			ok == true;
			rz = atof(argv[1]);
			argc--; argv++;
		}

		if((ok ==false) && (strcmp(argv[1], "-t") == 0)){
			argc--; argv++;
			ok == true;
			tx = atof(argv[1]);
			argc--; argv++;
			ty = atof(argv[1]);
			argc--; argv++;
			tz = atof(argv[1]);
			argc--; argv++;

		}
		if((ok ==false) && (strcmp(argv[1], "-origin") == 0)){
			argc--; argv++;
			ok == true;
			ox = atof(argv[1]);
			argc--; argv++;
			oy = atof(argv[1]);
			argc--; argv++;
			oz = atof(argv[1]);
			argc--; argv++;

		}

		if((ok ==false) && (strcmp(argv[1], "-sg") == 0)){
			argc--; argv++;
			ok == true;
			sg = atof(argv[1]);
			argc--; argv++;
		}
		if((ok ==false) && (strcmp(argv[1], "-inputVol") == 0)){
			argc--; argv++;
			ok == true;
			inputImageFile = argv[1];
			argc--; argv++;
		}
		if((ok ==false) && (strcmp(argv[1], "-logFileName") == 0)){
			argc--; argv++;
			ok == true;
			logfilename = argv[1];
			argc--; argv++;
		}

	}

	//Archivo de Registro LOG
	std::ofstream logregistro(logfilename);
		
	if(verbose)
	{
		if (inputImageFile)  
			logregistro  << "Input image: "  << inputImageFile  << std::endl;
			//std::cout << "Input image: "  << inputImageFile  << std::endl;
	}
	if(random_mode && semillaRandom){
		srand(time(NULL));	
	}

	//Imprimir el uso del comando
	menu();

	const     unsigned int   Dimension = 3;
	typedef   double  InputPixelType;
	typedef   unsigned short  OutputPixelType;

	typedef itk::Image< InputPixelType,  Dimension >   InputImageType;
	typedef itk::Image< OutputPixelType, Dimension >   OutputImageType;


	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	// Read the input image
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inputImageFile );
	reader->Update();

	
	for(int i=0; i<numberOfImages ; i++)
	{

		WriterType::Pointer writer = WriterType::New();
		
		//Resampler para el cambio del volumen
		typedef itk::ResampleImageFilter<InputImageType,OutputImageType> ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		
		//Tipo de Transformacion
		typedef itk::Similarity3DTransform< double >  SimilarityTransformType;
		SimilarityTransformType::Pointer similarityTransform = SimilarityTransformType::New();
		resample->SetTransform( similarityTransform );

		//Interpolador usado en el Resampler
		typedef itk::LinearInterpolateImageFunction<InputImageType, double >  InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		resample->SetInterpolator( interpolator );

		//Lectura de la imagen
		InputImageType::Pointer image = reader->GetOutput();


		//Image's spacing
		InputImageType::SpacingType spacing = image->GetSpacing();
		//Image's origin
		InputImageType::PointType origin;
		origin = image->GetOrigin();
		
		//In case of new origin
		if(ox != 0 && oy != 0 && oz != 0){
			InputImageType::PointType neworigin;
			neworigin[0] = 0.0;
			neworigin[1] = 0.0;
			neworigin[2] = 0.0;
			image->SetOrigin(neworigin);
		}	
		
		InputImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
		SimilarityTransformType::InputPointType center;
	
		//In case of new center of rotation	
		if(centerSpec){
			for(unsigned int j=0; j< Dimension; j++)
			{
				center[j] = reader->GetOutput()->GetOrigin()[j] + spacing[j]*size[j] / 2.0;
			}
			similarityTransform->SetCenter(center); 
		}
		//Initializing Similarity Transform
		similarityTransform->SetIdentity();

		//Parametros para el Versor 3D Rotation
		typedef SimilarityTransformType::VersorType VersorType;
		typedef VersorType::VectorType VectorType;
		VersorType rotation;
		VectorType axis;
		
		//Utilitario donde esta el conversor de rotacion euler a rotacion versor
		Utilitarios *util = new Utilitarios();
		//Procesamiento de parametros con modo Aleatorio o modo User
		
		if(random_mode){
			
			//Modo Aleatorio: parametros de transformacion definidos aleatoriamente
	
			
			//Versor 3D Rotation
			//rotaciones entre -10 y 10 grados sexag 
			//rx = rand()%(10 - (-10) + 1) +(-10);
			//ry = rand()%(10 - (-10) + 1) +(-10);
			//rz = rand()%(10 - (-10) + 1) +(-10);
			rx = ((double)rand()/RAND_MAX)*(10 - -10) + -10;
			ry = ((double)rand()/RAND_MAX)*(10 - -10) + -10;
			rz = ((double)rand()/RAND_MAX)*(10 - -10) + -10;

			util->convertEulerToVersor(rx,ry,rz,ax,ay,az,angle);
			
			//Warning no se considera el angulo en los parametros de rotacion		
			axis[0] = ax; axis[1] = ay; axis[2] = az;
			rotation.Set(axis,angle);
			similarityTransform->SetRotation(rotation);

			SimilarityTransformType::ParametersType similarityParameters;
			similarityParameters = similarityTransform->GetParameters();

			//Traslation Vector 
			//traslaciones entre -10 y 10 mm
			similarityParameters[3] = ((double)rand()/RAND_MAX)*(10 - -10) + -10;
			similarityParameters[4] = ((double)rand()/RAND_MAX)*(10 - -10) + -10;
			similarityParameters[5] = ((double)rand()/RAND_MAX)*(10 - -10) + -10;

			//Scale Factor (no podemos afectar mucho la escala)
			//sin escalas grandes obviarian informacion de la imagen
			//escala entre 1.0 y 1.5 
			similarityParameters[6] = ((double)rand()/RAND_MAX)*(1.5 -  1.0) + 1.0;
			//similarityParameters[6] = 1;

			similarityTransform->SetParameters(similarityParameters);
		
		}else{
			//Modo User: parametros de transformacion definidos por el usuario
			//Rotacion 3d de Similaridad
			util->convertEulerToVersor(rx,ry,rz,ax,ay,az,angle);	
			axis[0] = ax; axis[1] = ay; axis[2] = az;
			rotation.Set(axis,angle);
			similarityTransform->SetRotation(rotation);

			//Traslación 3D de Similaridad
			SimilarityTransformType::OutputVectorType translation;
			translation[0] = tx;
			translation[1] = ty;
			translation[2] = tz;
			similarityTransform->SetTranslation(translation);
			
			//Escala 3d de Similaridad
			SimilarityTransformType::ScaleType scale;
			scale = sg;
			similarityTransform->SetScale(scale);
		}



		// Initialize the resampler
		resample->SetSize(size);
		resample->SetOutputOrigin(origin);
		resample->SetOutputSpacing(spacing);
		resample->SetDefaultPixelValue(0);
		resample->SetOutputDirection( reader->GetOutput()->GetDirection());
		resample->SetInput(image);
		
		writer->SetInput( resample->GetOutput() );

		string fName = "";
		ostringstream fnameStream;
		fnameStream << i ;

		//Write the transform files
		itk::TransformFileWriterTemplate<double>::Pointer transformFileWriter =  itk::TransformFileWriterTemplate<double>::New();
		itksys::SystemTools::MakeDirectory( (fName + folderName + "/TransformFiles/").c_str() );
		string fileName = fName + folderName + "/TransformFiles/" + "transfSim_" +fnameStream.str() + ".txt";
		transformFileWriter->SetFileName(fileName.c_str());
		transformFileWriter->SetInput(similarityTransform);
		transformFileWriter->Update();

		itksys::SystemTools::MakeDirectory( (fName + folderName + "/Images/").c_str() );
		fName = fName + folderName + "/Images/" + "imagenDef_"+ fnameStream.str();
		fName += (Dimension == 2) ? ".png" : ".mha";

		writer->SetFileName( fName.c_str() );
		logregistro << "Writing " << fName.c_str() << std::endl; 
		//std::cout << "Writing " << fName.c_str() << std::endl;
		writer->Update();

		//Escribiendo Informacion General de Deformacion
		if(verbose){
			SimilarityTransformType::ParametersType finalParameters;
			finalParameters = similarityTransform->GetParameters();

			//La salida en Versor to Euler

			const double RotationAlongX = finalParameters[0]; // Convert radian to degree
			const double RotationAlongY = finalParameters[1];
			const double RotationAlongZ = finalParameters[2];
			const double TranslationAlongX = finalParameters[3];
			const double TranslationAlongY = finalParameters[4];
			const double TranslationAlongZ = finalParameters[5];
			const double EscalaXYZ = finalParameters[6];
			
			logregistro << "Deformacion Aplicada = " << std::endl;
			logregistro << " Rotation Degrees = [ " << rx << ", "<< ry << ", " << rz  << " ]"<<std::endl;
			logregistro << " Versor = [ " << RotationAlongX << ", "<< RotationAlongY << ", " << RotationAlongZ  << " ]"<<std::endl;
			logregistro << " Translation = [ " << TranslationAlongX  << ", "<<TranslationAlongY << ", "<< TranslationAlongZ << " ]"<<std::endl;
			logregistro << " Scale  = " << EscalaXYZ << std::endl;

			//std::cout << "Deformacion Aplicada = " << std::endl;
			//std::cout << " Rotation Degrees = [ " << rx << ", "<< ry << ", " << rz  << " ]"<<std::endl;
			//std::cout << " Versor = [ " << RotationAlongX << ", "<< RotationAlongY << ", " << RotationAlongZ  << " ]"<<std::endl;
			//std::cout << " Translation = [ " << TranslationAlongX  << ", "<<TranslationAlongY << ", "<< TranslationAlongZ << " ]"<<std::endl;
			//std::cout << " Scale  = " << EscalaXYZ << std::endl;



		}
		logregistro.close();

	}
	return EXIT_SUCCESS;
}

