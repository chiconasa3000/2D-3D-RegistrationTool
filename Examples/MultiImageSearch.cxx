
#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <itkCommand.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkExhaustiveOptimizer.h>
#include <itkTranslationTransform.h>

#include "itkGradientDifferenceMultiImageToImageMetric.h"
#include "itkNormalizedGradientCorrelationMultiImageToImageMetric.h"
#include "itkPatternIntensityMultiImageToImageMetric.h"
#include "itkMultiResolutionMultiImageToImageRegistrationMethod.h"
#include "itkPatchedRayCastInterpolateImageFunction.h"

#include <vector>
#include <sstream>


/** 
 * MultiImageSearch
 *  
 * This example takes a moving image and an arbitrary number of fixed images
 * and explores a metric's values for multiple parameters of translations 
 * along x, y and z. This program's results should be used to plot the 
 * metric's values and examine its shape. Transform parameters are always 
 * tested around the identity transform (i.e. no translation in any direction)
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
typedef RegistrationType* RegistrationPointer;
 
void Execute(itk::Object *caller, const itk::EventObject & event)
{
  RegistrationPointer registration = 
    dynamic_cast<RegistrationPointer>( caller );
  if( itk::IterationEvent().CheckEvent( & event ) )
    {
    registration->SetInitialTransformParametersOfNextLevel(
      registration->GetInitialTransformParameters() );
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
typedef itk::ExhaustiveOptimizer  OptimizerType;
typedef OptimizerType*            OptimizerPointer;
 
void Execute(itk::Object *caller, const itk::EventObject & event)
{
  OptimizerPointer optimizer = 
    dynamic_cast<OptimizerPointer>( caller );
  if( itk::IterationEvent().CheckEvent( & event ) )
    {
    std::cout << optimizer->GetCurrentPosition() << " " <<
      optimizer->GetCurrentValue() << std::endl;
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
  if( argc < 11 )
    {
    std::cerr << "Usage: " << argv[0] << " [movingImage] [N] [fixedImage1] "
      << "[focalPoint1_x] [focalPoint1_y] [focalPoint1_z] ... "
      << "[focalPointN_x] [focalPointN_y] [focalPointN_z] "
      << "[metricType: gd, gc, pi] <sigma> [gridRadius_x] [gridRadius_y] "
      << " [gridRadius_z] [gridStep]" << std::endl;
    std::cerr << "Parameter sigma must be only given if metricType is pi"
      << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimensions = 3;

  typedef itk::Image<short,Dimensions>   FixedImageType;
  typedef itk::Image<short,Dimensions>  MovingImageType;
  typedef itk::MultiResolutionMultiImageToImageRegistrationMethod< 
    FixedImageType,MovingImageType >    RegistrationType;

  RegistrationType::Pointer registration = RegistrationType::New();

  typedef RegistrationObserver<FixedImageType,
                       MovingImageType> RegistrationObserverType;
  
  RegistrationObserverType::Pointer registrationObserver = 
    RegistrationObserverType::New();
  registration->AddObserver( itk::IterationEvent(), registrationObserver );


//----------------------------------------------------------------------------
// Create the transform
//----------------------------------------------------------------------------
  typedef itk::TranslationTransform< double, Dimensions > TransformType;

  TransformType::Pointer transform = TransformType::New();
  registration->SetTransform( transform );


//----------------------------------------------------------------------------
// Load the moving image
//----------------------------------------------------------------------------
  typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    MovingImageType,
                                    MovingImageType  > MovingImagePyramidType;

  argc--;
  argv++;

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
  argc--;
  argv++;

  registration->SetMovingImage( movingReader->GetOutput() );

  MovingImagePyramidType::Pointer movingPyramidFilter = 
    MovingImagePyramidType::New();
  registration->SetMovingImagePyramid( movingPyramidFilter );


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
    fixedReader->SetFileName( argv[0] );
    try
      {
      fixedReader->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cout << e.GetDescription() << std::endl;
      return EXIT_FAILURE;
      }
    argc--;
    argv++;

    FixedImageConstPointer fixedImage = fixedReader->GetOutput();
    registration->AddFixedImage( fixedImage );

    registration->AddFixedImageRegion( fixedImage->GetBufferedRegion() );

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    FocalPointType focalPoint;
    focalPoint[0] = atof( argv[0] );
    focalPoint[1] = atof( argv[1] );
    focalPoint[2] = atof( argv[2] );
    argc -= 3;
    argv += 3;
    interpolator->SetFocalPoint( focalPoint );
    interpolator->SetTransform( transform );
    interpolator->SetThreshold( 0.0 );
    registration->AddInterpolator( interpolator );

    FixedImagePyramidType::Pointer fixedPyramidFilter = 
      FixedImagePyramidType::New();
    registration->AddFixedImagePyramid( fixedPyramidFilter );
    }


//----------------------------------------------------------------------------
// Create the multi metric
//----------------------------------------------------------------------------
  std::string multiMetricStr = argv[0];
  
  if ( multiMetricStr.compare("gd") == 0 )
    {
    itk::GradientDifferenceMultiImageToImageMetric<
      FixedImageType,MovingImageType>::Pointer gdMetric = 
      itk::GradientDifferenceMultiImageToImageMetric<
      FixedImageType,MovingImageType>::New();
    registration->SetMultiMetric( gdMetric );
    }
  else if(multiMetricStr.compare("gc") == 0 )
    {
    itk::NormalizedGradientCorrelationMultiImageToImageMetric<
      FixedImageType,MovingImageType>::Pointer gcMetric = 
      itk::NormalizedGradientCorrelationMultiImageToImageMetric<
      FixedImageType,MovingImageType>::New();
    registration->SetMultiMetric( gcMetric );
    }
  else if(multiMetricStr.compare("pi") == 0 )
    {
    itk::PatternIntensityMultiImageToImageMetric<
      FixedImageType,MovingImageType>::Pointer piMetric = 
      itk::PatternIntensityMultiImageToImageMetric<
      FixedImageType,MovingImageType>::New();
    
    // read the sigma argument too
    piMetric->SetSigma( atof(argv[1]) );
    argc--;
    argv++;

    registration->SetMultiMetric( piMetric );
    }
  else
    {
    std::cerr << "Incorrect multi metric chosen" << std::endl;
    std::cerr << "Allowed values are gd (gradient difference), gc (gradient"
      " correlation) and pi (pattern intensity)" << std::endl;
    return EXIT_FAILURE;
    }

  argc--;
  argv++;


//----------------------------------------------------------------------------
// Create the optimizer
//----------------------------------------------------------------------------  
  typedef itk::ExhaustiveOptimizer OptimizerType;

  OptimizerType::Pointer optimizer = OptimizerType::New();

  unsigned int ParTotal = transform->GetNumberOfParameters();
  
  OptimizerType::ScalesType scales( ParTotal );
  scales.Fill(itk::NumericTraits<OptimizerType::ScalesType::ValueType>::One);
  optimizer->SetScales( scales );

  OptimizerType::StepsType steps( ParTotal );
  for( unsigned int p=0; p<ParTotal; p++ )
    {
    steps[p] = atoi( argv[p] );
    }
  optimizer->SetNumberOfSteps( steps );
  argc -= 3;
  argv += 3;

  optimizer->SetStepLength( atof( argv[0] ) );
  argc--;
  argv++;

  OptimizerObserver::Pointer observer = OptimizerObserver::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  registration->SetOptimizer( optimizer );


//----------------------------------------------------------------------------
// Set the moving and fixed images' schedules
//----------------------------------------------------------------------------
  const unsigned int ResolutionLevels = 3;

  RegistrationType::ScheduleType fixedSchedule( ResolutionLevels,Dimensions );
  fixedSchedule[0][0] = 4;
  fixedSchedule[0][1] = 4;
  fixedSchedule[0][2] = 1;
  fixedSchedule[1][0] = 2;
  fixedSchedule[1][1] = 2;
  fixedSchedule[1][2] = 1;
  fixedSchedule[2][0] = 1;
  fixedSchedule[2][1] = 1;
  fixedSchedule[2][2] = 1;

  RegistrationType::ScheduleType movingSchedule( ResolutionLevels,Dimensions);
  movingSchedule[0][0] = 4;
  movingSchedule[0][1] = 4;
  movingSchedule[0][2] = 4;
  movingSchedule[1][0] = 2;
  movingSchedule[1][1] = 2;
  movingSchedule[1][2] = 2;
  movingSchedule[2][0] = 1;
  movingSchedule[2][1] = 1;
  movingSchedule[2][2] = 1;

  registration->SetSchedules( fixedSchedule, movingSchedule );
  

//----------------------------------------------------------------------------
// Start the registration
//----------------------------------------------------------------------------
  typedef RegistrationType::ParametersType ParametersType;

  ParametersType initialParameters( ParTotal );
  initialParameters.Fill(itk::NumericTraits<ParametersType::ValueType>::Zero);
  registration->SetInitialTransformParameters( initialParameters );

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
// End of the example
//----------------------------------------------------------------------------
  return EXIT_SUCCESS;
}
