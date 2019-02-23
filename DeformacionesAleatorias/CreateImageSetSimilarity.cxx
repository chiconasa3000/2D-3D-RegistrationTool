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

//#include "itkAffineTransform.h"
#include <itkSimilarity3DTransform.h>
#include "itkBSplineDeformableTransform.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkTransformFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"

#include <string>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include <itksys/SystemTools.hxx>

using namespace std;

int main( int argc, char * argv[] )
{
	if( argc < 4 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << "  inputImageFile  folderName numberOfImages <on/off (opt: for srand)>" << std::endl;
		return EXIT_FAILURE;
	}

	if(argc == 5)
	{
		/* initialize random seed: */
		srand(time(NULL));
	}
	const     unsigned int   Dimension = 3;
	typedef   double  InputPixelType;
	typedef   unsigned short  OutputPixelType;

	typedef itk::Image< InputPixelType,  Dimension >   InputImageType;
	typedef itk::Image< OutputPixelType, Dimension >   OutputImageType;


	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	// Read the input image
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );
	reader->Update();

	// Rescale the input image to 0-255
	/*typedef itk::RescaleIntensityImageFilter<InputImageType, InputImageType >  RescaleFilterType;
	RescaleFilterType::Pointer    rescaleFilter    = RescaleFilterType::New();
	rescaleFilter->SetInput(    reader->GetOutput() );
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
	rescaleFilter->Update();*/

	int numberOfImages = atoi(argv[3]);

	for(int i=0; i<numberOfImages ; i++)
	{

		WriterType::Pointer writer = WriterType::New();

		typedef itk::ResampleImageFilter<InputImageType,OutputImageType> ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();


		typedef itk::Similarity3DTransform< double >  SimilarityTransformType;
		SimilarityTransformType::Pointer similarityTransform = SimilarityTransformType::New();

		typedef itk::BSplineDeformableTransform< double,
			Dimension,
			3 >     BSplineTransformType;

		resample->SetTransform( similarityTransform );
		typedef itk::LinearInterpolateImageFunction<
			InputImageType, double >  InterpolatorType;

		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		resample->SetInterpolator( interpolator );


		// Set the parameters of the affine transform
		//Get the spacing
		InputImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
		//Get the origin
		//BSplineTransformType::OriginType origin;
		//origin = reader->GetOutput()->GetOrigin();

		InputImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

		/*SimilarityTransformType::FixedParametersType center; 

		  for(unsigned int j=0; j< Dimension; j++)
		  {
		  center[j] = origin[j] + spacing[j]*size[j] / 2.0;
		  }*/
		similarityTransform->SetIdentity();
		//Reemplazado por setfixedParameters que es donde setea el centro de rotacion
		//similarityTransform->SetCenter(center);
		SimilarityTransformType::ParametersType similarityParameters;
		similarityParameters = similarityTransform->GetParameters();

		//Considerando al formula de numeros aleatorios
		//x = ((double)rand()/RAND_MAX ) * (fmax - fmin) + fmin;
		//x = (rand()/100.0 ) (fmax) -(rand()/100.0)(fmin) + fmin;
		//rand()%100/100.0*0.05 - 0.5*0.05 + 1.0

		//Tres formas de disponer un nro aleatorio
		//(rand()%100/100.0 - 0.5)*0.3;
		//(rand()%100/100.0 - 0.5)*0.2;
		//1.0 + (rand()%100/100.0 - 0.5)*0.05;
		//(rand()%100/100.0 - 0.5)*15.0;

		const double dtr = (atan(1.0) * 4.0)/180.0;

		//Versor 3D Rotation
		//rotaciones entre -5 y -5 grados sexag 
		//similarityParameters[0] = dtr*(rand()%(5 - (-5) + 1) +(-5));
		//similarityParameters[1] = dtr*(rand()%(5 - (-5) + 1) +(-5));
		//similarityParameters[2] = dtr*(rand()%(5 - (-5) + 1) +(-5));
		similarityParameters[0] = dtr*(0);
		similarityParameters[1] = dtr*(0);
		similarityParameters[2] = dtr*(0);



		//Traslation Vector 
		//traslaciones entre -10 y 10 mm
		//similarityParameters[3] = rand()%(10 - (-10) + 1)+(-10);
		//similarityParameters[4] = rand()%(10 - (-10) + 1)+(-10);
		//similarityParameters[5] = rand()%(10 - (-10) + 1)+(-10);

		similarityParameters[3] = 0;
		similarityParameters[4] = 0;
		similarityParameters[5] = 0;


		//Scale Factor (no podemos afectar mucho la escala)
		//sin escalas grandes obviarian informacion de la imagen
		//escala entre 0.6 y 0.8 
		//similarityParameters[6] = ((double)rand()/RAND_MAX)*(0.8 -  0.6) + 0.6;
		similarityParameters[6] = 1;


		similarityTransform->SetParameters(similarityParameters);
		//similarityTransform->SetFixedParameters(center); 
		// Initialize the resampler
		// Get the size of the image
		size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

		resample->SetSize(size);
		resample->SetOutputOrigin(reader->GetOutput()->GetOrigin());
		resample->SetOutputSpacing(spacing);
		resample->SetDefaultPixelValue( 0 );
		resample->SetOutputDirection( reader->GetOutput()->GetDirection());
		resample->SetInput( reader->GetOutput() );
		writer->SetInput( resample->GetOutput() );

		
		string fname;
		ostringstream fnameStream;
		fnameStream << i ;


		//Write the transform files
		//itk::TransformFileWriter::Pointer  transformFileWriter = itk::TransformFileWriter::New();
		itk::TransformFileWriterTemplate<double>::Pointer transformFileWriter =  itk::TransformFileWriterTemplate<double>::New();
		itksys::SystemTools::MakeDirectory( (fname + argv[2] + "/TransformFiles/").c_str() );

		string fileName = fname + argv[2] + "/TransformFiles/" + "transfSim_" +fnameStream.str() + ".txt";
		transformFileWriter->SetFileName(fileName.c_str());
		//transformFileWriter->SetPrecision(12);
		transformFileWriter->SetInput(similarityTransform);
		transformFileWriter->Update();

		itksys::SystemTools::MakeDirectory( (fname+argv[2]+"/Images/").c_str() );

		fname = fname + argv[2] + "/Images/" + "imagenDef_"+ fnameStream.str();
		if(Dimension == 2)
		{
			fname += ".png";
		}
		else
		{
			fname += ".mha";
		}

		writer->SetFileName( fname.c_str() );
		std::cout << "Writing " << fname.c_str() << std::endl;
		writer->Update();

	}
	return EXIT_SUCCESS;
}

