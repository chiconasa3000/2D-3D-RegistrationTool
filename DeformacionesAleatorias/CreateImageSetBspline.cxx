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
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "BSplineDeformableTransformOpt.h"

#include "itkImageRegionIterator.h"
#include "itkNormalVariateGenerator.h"

#include "itkTransformFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"
    
#include <string>
#include <sstream>
#include <fstream>

#include <cstdlib>
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

    typedef itk::BSplineDeformableTransformOpt< double,
                                           Dimension,
                                           3 >     BSplineTransformType;

    BSplineTransformType::Pointer bsplineTransform = BSplineTransformType::New();


    typedef BSplineTransformType::RegionType RegionType;
    RegionType bsplineRegion;
    RegionType::SizeType   gridSizeOnImage;
    RegionType::SizeType   gridBorderSize;
    RegionType::SizeType   totalGridSize;

    gridSizeOnImage.Fill( 8 );
    gridBorderSize.Fill( 3 );    // Border for spline order = 3 ( 1 lower, 2 upper )
    totalGridSize = gridSizeOnImage + gridBorderSize;

    bsplineRegion.SetSize( totalGridSize );

    typedef BSplineTransformType::SpacingType SpacingType;
    SpacingType bsplineSpacing = reader->GetOutput()->GetSpacing();

    typedef BSplineTransformType::OriginType OriginType;
    OriginType bsplineOrigin = reader->GetOutput()->GetOrigin();;

    InputImageType::SizeType ImageSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

    for(unsigned int r=0; r<Dimension; r++)
    {
      bsplineSpacing[r] *= floor( static_cast<double>(ImageSize[r] )  /
          static_cast<double>(gridSizeOnImage[r] ) );
      bsplineOrigin[r]  -=  bsplineSpacing[r];
    }

    bsplineTransform->SetGridSpacing( bsplineSpacing );
    bsplineTransform->SetGridOrigin( bsplineOrigin );
    bsplineTransform->SetGridRegion( bsplineRegion );
  

    typedef BSplineTransformType::ParametersType     ParametersType;

    const unsigned int numberOfParameters =
                     bsplineTransform->GetNumberOfParameters();
  
    ParametersType bsplineParameters( numberOfParameters );

    bsplineParameters.Fill( 0.0 );

    //Randomly set the parameters
    for(unsigned int j=0; j<bsplineParameters.GetSize(); j++)
    {
      bsplineParameters[j] = (rand()%203/203.0 - 0.5)*23.0;;
    }


    bsplineTransform->SetParameters( bsplineParameters );

    resample->SetTransform( bsplineTransform );
    typedef itk::LinearInterpolateImageFunction<
                             InputImageType, double >  InterpolatorType;

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resample->SetInterpolator( interpolator );

    // Initialize the resampler
    // Get the size of the image
    InputImageType::SizeType   size;
    size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();

    //Get the spacing
    InputImageType::SpacingType spacing;
    spacing = reader->GetOutput()->GetSpacing();
    //Get the origin
    BSplineTransformType::OriginType origin;
    origin = reader->GetOutput()->GetOrigin();
    resample->SetSize(size);
    resample->SetOutputOrigin(origin);
    resample->SetOutputSpacing(spacing);
    resample->SetOutputDirection( reader->GetOutput()->GetDirection());
    resample->SetDefaultPixelValue( 0 );
      
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
    transformFileWriter->SetInput(bsplineTransform);
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

