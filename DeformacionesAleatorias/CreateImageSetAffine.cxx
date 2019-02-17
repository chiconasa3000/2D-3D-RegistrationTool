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

#include "itkAffineTransform.h"
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
  typedef itk::RescaleIntensityImageFilter<InputImageType, InputImageType >  RescaleFilterType;
  RescaleFilterType::Pointer    rescaleFilter    = RescaleFilterType::New();
  rescaleFilter->SetInput(    reader->GetOutput() );
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  int numberOfImages = atoi(argv[3]);

  for(int i=0; i<numberOfImages ; i++)
  {

    WriterType::Pointer writer = WriterType::New();
   
    typedef itk::ResampleImageFilter<InputImageType,OutputImageType> ResampleFilterType;
    ResampleFilterType::Pointer resample = ResampleFilterType::New();

    typedef itk::AffineTransform< double, Dimension >  AffineTransformType;
    AffineTransformType::Pointer affineTransform = AffineTransformType::New();

    typedef itk::BSplineDeformableTransform< double,
                                           Dimension,
                                           3 >     BSplineTransformType;

    resample->SetTransform( affineTransform );
    typedef itk::LinearInterpolateImageFunction<
                             InputImageType, double >  InterpolatorType;

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resample->SetInterpolator( interpolator );
    

    // Set the parameters of the affine transform
    //Get the spacing
    InputImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
    //Get the origin
    BSplineTransformType::OriginType origin;
    origin = reader->GetOutput()->GetOrigin();
    
    InputImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

    AffineTransformType::InputPointType center;
    for(unsigned int j=0; j< Dimension; j++)
    {
      center[j] = origin[j] + spacing[j]*size[j] / 2.0;
    }
    affineTransform->SetIdentity();
    affineTransform->SetCenter(center);
    AffineTransformType::ParametersType affineParameters;
    affineParameters = affineTransform->GetParameters();

    //Considerando al formula de numeros aleatorios
    //x = ((double)rand()/RAND_MAX ) * (fmax - fmin) + fmin;
    //x = (rand()/100.0 ) (fmax) -(rand()/100.0)(fmin) + fmin;
    //rand()%100/100.0*0.05 - 0.5*0.05 + 1.0

    affineParameters[0] = 1.0 + (rand()%100/100.0 - 0.5)*0.05;
    affineParameters[1] = (rand()%100/100.0 - 0.5)*0.3;
    affineParameters[2] = (rand()%100/100.0 - 0.5)*0.2;
      
    affineParameters[3] = (rand()%100/100.0 - 0.5)*0.3;
    affineParameters[4] = 1.0 + (rand()%100/100.0 - 0.5)*0.05;
    affineParameters[5] = (rand()%100/100.0 - 0.5)*0.3;
      
    affineParameters[6] = (rand()%100/100.0 - 0.5)*0.2;
    affineParameters[7] = (rand()%100/100.0 - 0.5)*0.3;
    affineParameters[8] = 1.0 + (rand()%100/100.0 - 0.5)*0.05;
      
    affineParameters[9] = (rand()%100/100.0 - 0.5)*15.0;
    affineParameters[10] = (rand()%100/100.0 - 0.5)*15.0;
    affineParameters[11] = (rand()%100/100.0 - 0.5)*15.0;
    
    affineTransform->SetParameters(affineParameters);
      
    // Initialize the resampler
    // Get the size of the image
    size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

    resample->SetSize(size);
    resample->SetOutputOrigin(origin);
    resample->SetOutputSpacing(spacing);
    resample->SetDefaultPixelValue( 0 );
    resample->SetOutputDirection( reader->GetOutput()->GetDirection());
    resample->SetInput( rescaleFilter->GetOutput() );
    writer->SetInput( resample->GetOutput() );
      

    string fname;
    ostringstream fnameStream;
    fnameStream << i ;


    //Write the transform files
    //itk::TransformFileWriter::Pointer  transformFileWriter = itk::TransformFileWriter::New();
    itk::TransformFileWriterTemplate<double>::Pointer transformFileWriter =  itk::TransformFileWriterTemplate<double>::New();
    itksys::SystemTools::MakeDirectory( (fname + argv[2] + "/TransformFiles/").c_str() );

    string fileName = fname + argv[2] + "/TransformFiles/" + fnameStream.str() + ".txt";
    transformFileWriter->SetFileName(fileName.c_str());
    //transformFileWriter->SetPrecision(12);
    transformFileWriter->SetInput(affineTransform);
    transformFileWriter->Update();

    itksys::SystemTools::MakeDirectory( (fname+argv[2]+"/Images/").c_str() );
    
    fname = fname + argv[2] + "/Images/" + fnameStream.str();
    if(Dimension == 2)
    {
      fname += ".png";
    }
    else
    {
      fname += ".mhd";
    }
    
    writer->SetFileName( fname.c_str() );
    std::cout << "Writing " << fname.c_str() << std::endl;
    writer->Update();

  }
  return EXIT_SUCCESS;
}

