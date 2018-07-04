
#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <itkCommand.h>
//#include <itkEuler3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkFRPROptimizer.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

#include "itkNormalizedGradientCorrelationMultiImageToImageMetric.h"
#include "itkMultiResolutionMultiImageToImageRegistrationMethod.h"
#include "itkPatchedRayCastInterpolateImageFunction.h"

#include <iostream>


/** 
 * MultiImageRegistration
 *  
 * This example is a basic application which takes a moving 3D image and 
 * registers it to a given number of fixed 2D projections. Registration is 
 * done maximising the Normalized Gradient Correlation metric with a rigid 
 * transform. Three resolution levels are used, using downsampling factors of 
 * 4 and 2 for the first and second levels and the original resolutions for 
 * the last one.
 * 
 */ 


//----------------------------------------------------------------------------
// Registration observer class
//----------------------------------------------------------------------------
template <class TFixedImage, class TMovingImage>
class RegistrationObserver : public itk::Command 
{
public:
typedef  RegistrationObserver<TFixedImage,TMovingImage> Self;

typedef  itk::Command             Superclass;
typedef  itk::SmartPointer<Self>  Pointer;
itkNewMacro( Self );

protected:
RegistrationObserver() {};

public:
typedef itk::MultiResolutionMultiImageToImageRegistrationMethod<
                                               TFixedImage,
                                               TMovingImage> RegistrationType;
typedef RegistrationType*   RegistrationPointer;
typedef itk::FRPROptimizer  OptimizerType;
 
void Execute(itk::Object *caller, const itk::EventObject & event)
{
  RegistrationPointer registration = 
    dynamic_cast<RegistrationPointer>( caller );
  if( itk::IterationEvent().CheckEvent( & event ) )
    {
    std::cout << std::endl << "Resolution level " << 
      registration->GetCurrentLevel() << std::endl;
//Half step length and tolerance in each level
    OptimizerType* optimizer = dynamic_cast<OptimizerType*>(
      registration->GetOptimizer() );

    optimizer->SetStepLength( 0.5*optimizer->GetStepLength() );
    std::cout << "StepLength set to " << optimizer->GetStepLength() 
      << std::endl;
    optimizer->SetStepTolerance( 0.5*optimizer->GetStepTolerance() );
    std::cout << "StepTolerance set to " << optimizer->GetStepTolerance() 
      << std::endl;
    }
}

void Execute(const itk::Object* caller, const itk::EventObject & event)
{}

};


//----------------------------------------------------------------------------
// Optimizer observer class
//----------------------------------------------------------------------------
class OptimizerObserver : public itk::Command 
{
public:
typedef  OptimizerObserver     Self;
typedef  itk::Command             Superclass;
typedef  itk::SmartPointer<Self>  Pointer;
itkNewMacro( Self );

protected:
OptimizerObserver() {};

public:
typedef itk::FRPROptimizer  OptimizerType;
typedef OptimizerType*      OptimizerPointer;
 
void Execute(itk::Object *caller, const itk::EventObject & event)
{
  OptimizerPointer optimizer = 
    dynamic_cast<OptimizerPointer>( caller );
  if( itk::IterationEvent().CheckEvent( & event ) )
    {
    std::cout << "Iteration " << optimizer->GetCurrentIteration()
      << "/" << optimizer->GetMaximumIteration() << " Position: " << 
      optimizer->GetCurrentPosition() << " Value: " << 
      optimizer->GetCurrentCost() << std::endl;
    }
  else if( itk::StartEvent().CheckEvent( & event ) )
    {
    std::cout << "Optimization started ..." << std::endl;
    }
  else if( itk::EndEvent().CheckEvent( & event ) )
    {
    std::cout << "Optimization ended." << std::endl;
    }
}


void Execute(const itk::Object* caller, const itk::EventObject & event)
{}

};


//----------------------------------------------------------------------------
// Main function
//----------------------------------------------------------------------------
int main(int argc, char* argv[] )
{
  if( argc < 13 )
    {
    std::cerr << "Usage: " << argv[0] << " [movingImage] [N] [fixedImage1] "
      << "[focalPoint1_x] [focalPoint1_y] [focalPoint1_z] ... "
      << "[focalPointN_x] [focalPointN_y] [focalPointN_z] "
      << "[tolerance] [stepsize]"
      << "[levels_n] [schedule_1,schedule_2...,schedule_n]"
      << "[outputDirectory] <inputTransform>" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimensions = 3;

  typedef itk::Image<float,Dimensions>  FixedImageType;
  typedef itk::Image<float,Dimensions>  MovingImageType;
  typedef itk::MultiResolutionMultiImageToImageRegistrationMethod< 
    FixedImageType,MovingImageType >    RegistrationType;

  RegistrationType::Pointer registration = RegistrationType::New();

  typedef RegistrationObserver<FixedImageType,
                       MovingImageType> RegistrationObserverType;
  
  RegistrationObserverType::Pointer registrationObserver = 
    RegistrationObserverType::New();
  registration->AddObserver( itk::IterationEvent(), registrationObserver );

  argc--;
  argv++;


//----------------------------------------------------------------------------
// Create the transform
//----------------------------------------------------------------------------
  typedef itk::Similarity3DTransform< double > TransformType;

  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  registration->SetTransform( transform );


//----------------------------------------------------------------------------
// Load the moving image
//----------------------------------------------------------------------------
  typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    MovingImageType,
                                    MovingImageType  > MovingImagePyramidType;

  MovingImageReaderType::Pointer movingReader = MovingImageReaderType::New();
  movingReader->SetFileName( argv[0] );
  try
    {
    movingReader->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  registration->SetMovingImage( movingReader->GetOutput() );

  MovingImagePyramidType::Pointer movingPyramidFilter = 
    MovingImagePyramidType::New();
  registration->SetMovingImagePyramid( movingPyramidFilter );

  argc--;
  argv++;


//----------------------------------------------------------------------------
// Load the fixed images
//----------------------------------------------------------------------------
  typedef FixedImageType::ConstPointer          FixedImageConstPointer;
  typedef itk::ImageFileReader<FixedImageType>  FixedImageReaderType;

  typedef itk::PatchedRayCastInterpolateImageFunction<MovingImageType,
                                                     double> InterpolatorType;
  typedef InterpolatorType::InputPointType FocalPointType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    FixedImageType,
                                    FixedImageType  >  FixedImagePyramidType;

  const unsigned int FImgTotal = atoi( argv[0] );
  argc--;
  argv++;

  for( unsigned int f=0; f<FImgTotal; f++ )
    {
    FixedImageReaderType::Pointer fixedReader = FixedImageReaderType::New();
    fixedReader->SetFileName( argv[ 0 ] );
    try
      {
      fixedReader->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cout << e.GetDescription() << std::endl;
      return EXIT_FAILURE;
      }

    FixedImageConstPointer fixedImage = fixedReader->GetOutput();
    registration->AddFixedImage( fixedImage );

    registration->AddFixedImageRegion( fixedImage->GetBufferedRegion() );

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    FocalPointType focalPoint;
    focalPoint[0] = atof( argv[ 1 ] );
    focalPoint[1] = atof( argv[ 2 ] );
    focalPoint[2] = atof( argv[ 3 ] );
    interpolator->SetFocalPoint( focalPoint );
    interpolator->SetTransform( transform );
    interpolator->SetThreshold( 0.0 );
    registration->AddInterpolator( interpolator );

    FixedImagePyramidType::Pointer fixedPyramidFilter = 
      FixedImagePyramidType::New();
    registration->AddFixedImagePyramid( fixedPyramidFilter );

    argc -= 4;
    argv += 4;
    }


//----------------------------------------------------------------------------
// Create the multi metric
//----------------------------------------------------------------------------
  typedef itk::NormalizedGradientCorrelationMultiImageToImageMetric<
    FixedImageType, MovingImageType> MultiMetricType;

  MultiMetricType::Pointer multiMetric = MultiMetricType::New();
  registration->SetMultiMetric( multiMetric );


//----------------------------------------------------------------------------
// Create the optimizer
//----------------------------------------------------------------------------  
  typedef itk::FRPROptimizer OptimizerType;

  OptimizerType::Pointer optimizer = OptimizerType::New();

  unsigned int ParTotal = transform->GetNumberOfParameters();
  
  OptimizerType::ScalesType scales( ParTotal );
  scales.Fill(itk::NumericTraits<OptimizerType::ScalesType::ValueType>::One);
  scales[0] = 25.0;
  scales[1] = 25.0;
  scales[2] = 25.0;
  optimizer->SetScales( scales );

  optimizer->SetMaximize( true );
  optimizer->SetMaximumIteration( 100 );
  optimizer->SetMaximumLineIteration( 10 );
  optimizer->SetValueTolerance( 1e-3 );
  optimizer->SetUseUnitLengthGradient( true );
  optimizer->SetToPolakRibiere();
  optimizer->SetCatchGetValueException( true );
  optimizer->SetMetricWorstPossibleValue( 
    -itk::NumericTraits<MultiMetricType::MeasureType>::infinity() );

//These parameters will be halved at the begginng of each resolution level
  optimizer->SetStepTolerance( atof(argv[0]) );
  optimizer->SetStepLength( atof(argv[1]) );
  argc-=2;
  argv+=2;

  OptimizerObserver::Pointer optimizerObserver = OptimizerObserver::New();
  optimizer->AddObserver( itk::IterationEvent(), optimizerObserver );
  optimizer->AddObserver( itk::StartEvent(), optimizerObserver );
  optimizer->AddObserver( itk::EndEvent(), optimizerObserver );

  registration->SetOptimizer( optimizer );


//----------------------------------------------------------------------------
// Set the moving and fixed images' schedules
//----------------------------------------------------------------------------
  const unsigned int ResolutionLevels =atoi(argv[0]);
  argc--;
  argv++;

  RegistrationType::ScheduleType movingSchedule( ResolutionLevels,Dimensions); 
  RegistrationType::ScheduleType fixedSchedule( ResolutionLevels,Dimensions );  
  
  int schedulesScales[ResolutionLevels];
  for(int i=0; i<ResolutionLevels; i++){
  	schedulesScales[i] = atoi(argv[i]);
  }
  argc-= ResolutionLevels;
  argv+= ResolutionLevels;
  
  std::cout<<"::Schedule::"<<std::endl;  

  for(int i=0; i<ResolutionLevels; i++){
	std::cout<<"Nivel: "<<i<<std::endl;
	for(int j=0; j<Dimensions; j++){
		if(j!=2){
			fixedSchedule[i][j]=schedulesScales[i];
			movingSchedule[i][j]=schedulesScales[i];
		}else{
			fixedSchedule[i][j]=1;
			movingSchedule[i][j]=schedulesScales[i];
		}
		
		std::cout<<"[x,y,z]: "<<"["<<fixedSchedule[i][j]<<", ";
	}
	std::cout<<"]";		
  }
 
  registration->SetSchedules( fixedSchedule, movingSchedule );


//----------------------------------------------------------------------------
// Get the output directory
//----------------------------------------------------------------------------
  std::string outDir = argv[0];

  argc--;
  argv++;


//----------------------------------------------------------------------------
// Read the transform file, if given
//----------------------------------------------------------------------------

  if( argc > 0 ) // if there are still arguments in the command line
    {
    itk::TransformFileReader::Pointer transformReader = 
      itk::TransformFileReader::New();
    transformReader->SetFileName( argv[0] );
    argc--;
    argv++;

    std::cout << "Attempting to read transform file " << 
      transformReader->GetFileName() << std::endl;
    try
      {
      transformReader->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cout << e.GetDescription() << std::endl;
      return EXIT_FAILURE;
      }

    itk::TransformFileReader::TransformListType*
      transformList = transformReader->GetTransformList();
    
    if( transformList->size() > 1 )
      {
      std::cout << "Only one transform expected from file. "
        << "Using first available transform" << std::endl;
      }
    
    itk::TransformFileReader::TransformPointer baseTransform =
      transformList->front();
    std::string className = baseTransform->GetNameOfClass();
    if( className.compare("Similarity3DTransform") != 0 )
      {
      std::cerr << "Transform class must be Similarity3DTransform." << std::endl;
      std::cerr << "Found " << className << " instead."
        << std::endl;
      return EXIT_FAILURE;
      }
    transform->SetFixedParameters( baseTransform->GetFixedParameters() );
    transform->SetParameters( baseTransform->GetParameters() );
    }

  registration->SetInitialTransformParameters( transform->GetParameters() );


//----------------------------------------------------------------------------
// Save the initial matrix file
//----------------------------------------------------------------------------
  std::string inputMatrixPath = outDir;
  inputMatrixPath.append( "/inMatrix.txt" );

  std::ofstream inputMatrixFile;
  std::cout << "Saving input matrix file " << inputMatrixPath << std::endl;
  inputMatrixFile.open( inputMatrixPath.c_str() );
  if(  inputMatrixFile.is_open() )
    {
    TransformType::MatrixType m = transform->GetMatrix();
    TransformType::OffsetType off = transform->GetOffset();
    inputMatrixFile << m[0][0] << " " << m[0][1] << " " << m[0][2] << " "
      << off[0] << std::endl;
    inputMatrixFile << m[1][0] << " " << m[1][1] << " " << m[1][2] << " " 
      << off[1] << std::endl;
    inputMatrixFile << m[2][0] << " " << m[2][1] << " " << m[2][2] << " " 
      << off[2] << std::endl;
    inputMatrixFile << "0.0 0.0 0.0 1.0" << std::endl;
    inputMatrixFile.close();
    }
  else
    {
    std::cerr << "Could not open input matrix file " << inputMatrixPath << std::endl;
    }


//----------------------------------------------------------------------------
// Start the registration
//----------------------------------------------------------------------------

  try
    {
    registration->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }


//----------------------------------------------------------------------------
// Save the transform files
//----------------------------------------------------------------------------

  std::string outputTransformPath = outDir;
  outputTransformPath.append( "/outTransform.txt" );

  itk::TransformFileWriter::Pointer transformWriter = 
    itk::TransformFileWriter::New();
  transformWriter->SetInput( transform );
  transformWriter->SetFileName( outputTransformPath.c_str() );

  std::cout << "Saving output transform file " << transformWriter->GetFileName() 
    << std::endl;
  try
    {
    transformWriter->Update();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cout << e.GetDescription() << std::endl;
    }

  std::string outputMatrixPath = outDir;
  outputMatrixPath.append( "/outMatrix.txt" );

  std::ofstream outputMatrixFile;
  std::cout << "Saving output matrix file " << outputMatrixPath << std::endl;
  outputMatrixFile.open( outputMatrixPath.c_str() );
  if(  outputMatrixFile.is_open() )
    {
    TransformType::MatrixType m = transform->GetMatrix();
    TransformType::OffsetType off = transform->GetOffset();
    outputMatrixFile << m[0][0] << " " << m[0][1] << " " << m[0][2] << " "
      << off[0] << std::endl;
    outputMatrixFile << m[1][0] << " " << m[1][1] << " " << m[1][2] << " " 
      << off[1] << std::endl;
    outputMatrixFile << m[2][0] << " " << m[2][1] << " " << m[2][2] << " " 
      << off[2] << std::endl;
    outputMatrixFile << "0.0 0.0 0.0 1.0" << std::endl;
    outputMatrixFile.close();
    }
  else
    {
    std::cerr << "Could not open output matrix file " << outputMatrixPath << std::endl;
    }


//----------------------------------------------------------------------------
// Save projected images
//----------------------------------------------------------------------------
  typedef itk::ResampleImageFilter<MovingImageType,
                                   FixedImageType>  ResamplerType;
  typedef itk::SubtractImageFilter<FixedImageType,FixedImageType,
                                   FixedImageType>  SubtracterType;
  typedef itk::ImageFileWriter<FixedImageType>      WriterType;


  for( unsigned int f=0; f<FImgTotal; f++ )
    {
    ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetInput( registration->GetMovingImage() );
    resampler->SetInterpolator( registration->GetMultiInterpolator()[f] );
    resampler->SetTransform( transform );
    resampler->SetUseReferenceImage( true );
    resampler->SetReferenceImage( registration->GetFixedMultiImage()[f] );

    typedef itk::RescaleIntensityImageFilter< MovingImageType, MovingImageType > RescaleFilterType;

    RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();
    rescaler1->SetOutputMinimum(   0 );
    rescaler1->SetOutputMaximum( 255 );
    rescaler1->SetInput( resampler->GetOutput());
    rescaler1->Update();

    std::stringstream strStream;
    strStream << outDir << "/projection";
    strStream.width(2);
    strStream.fill('0');
    strStream << f;
    strStream.width(0);
    strStream << ".mha";

    WriterType::Pointer projectionWriter = WriterType::New();
    projectionWriter->SetFileName( strStream.str().c_str() );
    //projectionWriter->SetInput( resampler->GetOutput() );
    projectionWriter->SetInput( rescaler1->GetOutput() );

    std::cout << "Attempting to write projection file " << 
      projectionWriter->GetFileName() << std::endl;
    try
      {
      projectionWriter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << e.GetDescription() << std::endl;
      }

    SubtracterType::Pointer subtracter = SubtracterType::New();
    subtracter->SetInput1( registration->GetFixedMultiImage()[f] );
    //subtracter->SetInput2( resampler->GetOutput() );
    subtracter->SetInput2( rescaler1->GetOutput() );
    
    std::stringstream strStream2;
    strStream2 << outDir << "/subtraction";
    strStream2.width(2);
    strStream2.fill('0');
    strStream2 << f;
    strStream2.width(0);
    strStream2 << ".mha";

    WriterType::Pointer subtractionWriter = WriterType::New();
    subtractionWriter->SetFileName( strStream2.str().c_str() );
    subtractionWriter->SetInput( subtracter->GetOutput() );

    std::cout << "Attempting to write subtraction file " << 
      subtractionWriter->GetFileName() << std::endl;
    try
      {
      subtractionWriter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << e.GetDescription() << std::endl;
      }

    }


//----------------------------------------------------------------------------
// End of the example
//----------------------------------------------------------------------------
  return EXIT_SUCCESS;
}

