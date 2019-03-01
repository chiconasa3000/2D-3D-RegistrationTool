#ifndef __itkMultiResolutionMultiImageToImageRegistrationMethod_txx
#define __itkMultiResolutionMultiImageToImageRegistrationMethod_txx

#include "itkMultiResolutionMultiImageToImageRegistrationMethod.h"
#include <itkRecursiveMultiResolutionPyramidImageFilter.h>
#include <sstream>
#include <string.h>
#include <fstream>
#include <iostream>
namespace itk
{

/**
 * Constructor
 */
template < typename TFixedImage, typename TMovingImage >
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::MultiResolutionMultiImageToImageRegistrationMethod()
{
    m_NumberOfLevels = 1;
    m_CurrentLevel = 0;

    m_Stop = false;

    m_ScheduleSpecified = false;
    m_NumberOfLevelsSpecified = false;

    m_InitialTransformParametersOfNextLevel = ParametersType(1);

    m_InitialTransformParametersOfNextLevel.Fill( 0.0f );
}


/*
 * Initialize by setting the interconnects between components.
 */
template < typename TFixedImage, typename TMovingImage >
void
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::Initialize() throw (ExceptionObject)
{
    Superclass::Initialize();

    // Setup the metric
    this->m_MultiMetric->SetMovingImage( m_MovingImagePyramid->GetOutput( m_CurrentLevel ) );

    const unsigned int FImgTotal = this->GetNumberOfFixedImages();

    // Create the current vector of fixed images
    while( this->m_CurrentFixedMultiImage.size() > 0 )
    {
        this->m_CurrentFixedMultiImage.pop_back();
    }
    /*
    for(unsigned int f=0; f<FImgTotal; f++ )
    {
        this->m_CurrentFixedMultiImage.push_back( this->m_FixedMultiImagePyramid[f]->GetOutput( m_CurrentLevel ) );

        //----------------------------------------------------------------------
        //save the moving image projection with the metric has already applied
        //---------------------------------------------------------------------

        //the moving image is updating with the last transformation
        typename ResamplerType::Pointer resampler = ResamplerType::New();
        resampler->SetInput( m_MovingImagePyramid->GetOutput( m_CurrentLevel ) );
        resampler->SetInterpolator(this->m_MultiInterpolator[f] );
        resampler->SetTransform( this->m_Transform );
        resampler->SetUseReferenceImage( true );
        resampler->SetReferenceImage(this->m_FixedMultiImagePyramid[f]->GetOutput( m_CurrentLevel ) );


        //threshold the image with the correct intensity values
        typename RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();
        rescaler1->SetOutputMinimum(   0 );
        rescaler1->SetOutputMaximum( 255 );
        rescaler1->SetInput( resampler->GetOutput());
        rescaler1->Update();

        std::stringstream strStream;
        strStream <<"current_projection";
        strStream.width(1);
        strStream << m_CurrentLevel;
        strStream << f;
        strStream.width(0);
        strStream << ".mha";

        //the projection is condition
        typename WriterType::Pointer projectionWriter = WriterType::New();
        projectionWriter->SetFileName( strStream.str().c_str() );
        //projectionWriter->SetInput( resampler->GetOutput() );
        projectionWriter->SetInput( resampler->GetOutput() );

        std::cout << "Attempting to write projection file " << projectionWriter->GetFileName() << std::endl;
        try
        {
            projectionWriter->Update();
        }
        catch( itk::ExceptionObject & e )
        {
            std::cerr << e.GetDescription() << std::endl;
        }

        //the difference is between the projection and the current fixed image in this case the last level
        typename SubtracterType::Pointer subtracter = SubtracterType::New();
        subtracter->SetInput1( this->m_FixedMultiImagePyramid[f]->GetOutput( m_CurrentLevel ) );
        //subtracter->SetInput2( resampler->GetOutput() );
        subtracter->SetInput2( rescaler1->GetOutput() );

        std::stringstream strStream2;
        strStream2 <<"current_subtraction";
        strStream2.width(1);
        strStream2 << m_CurrentLevel;
        strStream2 << f;
        strStream2.width(0);
        strStream2 << ".mha";

        typename WriterType::Pointer subtractionWriter = WriterType::New();
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
    }*/

    // Create the current vector of fixed images' regions
    while( this->m_CurrentFixedMultiImageRegion.size() > 0 )
    {
        this->m_CurrentFixedMultiImageRegion.pop_back();
    }
    for(unsigned int f=0; f<FImgTotal; f++ )
    {
        this->m_CurrentFixedMultiImageRegion.push_back(this->m_FixedMultiImageRegionPyramid[f][m_CurrentLevel] );
    }

    this->m_MultiMetric->SetFixedMultiImage( this->m_CurrentFixedMultiImage );
    this->m_MultiMetric->SetTransform( this->m_Transform );
    this->m_MultiMetric->SetMultiInterpolator( this->m_MultiInterpolator );
    this->m_MultiMetric->SetFixedMultiImageRegion( this->m_CurrentFixedMultiImageRegion );
    this->m_MultiMetric->Initialize();

    // Setup the optimizer
    this->m_Optimizer->SetCostFunction( this->m_MultiMetric );
    this->m_Optimizer->SetInitialPosition( m_InitialTransformParametersOfNextLevel );

    // Connect the transform to the Decorator.
    TransformOutputType * transformOutput =
            static_cast< TransformOutputType * >( this->ProcessObject::GetOutput(0) );

    transformOutput->Set( this->m_Transform.GetPointer() );

}


/*
 * Stop the Registration Process
 */
template < typename TFixedImage, typename TMovingImage >
void
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::StopRegistration( void )
{
    m_Stop = true;
}

/**
 * Set the schedules for the fixed and moving image pyramid
 */
template < typename TFixedImage, typename TMovingImage >
void
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::SetSchedules( const ScheduleType & fixedImagePyramidSchedule,
                const ScheduleType & movingImagePyramidSchedule )
{
    if( m_NumberOfLevelsSpecified )
    {
        itkExceptionMacro( "SetSchedules should not be used "
                           << "if numberOfLevelves are specified using SetNumberOfLevels" );
    }
    m_FixedImagePyramidSchedule = fixedImagePyramidSchedule;
    m_MovingImagePyramidSchedule = movingImagePyramidSchedule;
    m_ScheduleSpecified = true;

    //Set the number of levels based on the pyramid schedule specified
    if ( m_FixedImagePyramidSchedule.rows() !=
         m_MovingImagePyramidSchedule.rows())
    {
        itkExceptionMacro("The specified schedules contain unequal number of levels");
    }
    else
    {
        m_NumberOfLevels = m_FixedImagePyramidSchedule.rows();
    }

    this->Modified();
}

/**
 * Set the number of levels
 */
template < typename TFixedImage, typename TMovingImage >
void
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::SetNumberOfLevels( unsigned long numberOfLevels )
{
    if( m_ScheduleSpecified )
    {
        itkExceptionMacro( "SetNumberOfLevels should not be used "
                           << "if schedules have been specified using SetSchedules method " );
    }

    m_NumberOfLevels = numberOfLevels;
    m_NumberOfLevelsSpecified = true;
    this->Modified();
}

/**
 * Prepare the image pyramids.
 */
template < typename TFixedImage, typename TMovingImage >
void
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::PreparePyramids( void )
{
    if( !this->m_Transform )
    {
        itkExceptionMacro(<<"Transform is not present");
    }

    m_InitialTransformParametersOfNextLevel = this->m_InitialTransformParameters;

    if ( m_InitialTransformParametersOfNextLevel.Size() != this->m_Transform->GetNumberOfParameters() )
    {
        itkExceptionMacro(<<"Size mismatch between initial parameter and transform");
    }

    // Sanity checks
    const unsigned int FImgTotal = this->GetNumberOfFixedImages();
    if( FImgTotal == 0 )
    {
        itkExceptionMacro(<<"Fixed multi image is empty");
    }

    if( !this->m_MovingImage )
    {
        itkExceptionMacro(<<"Moving image is not present");
    }

    // Check the pyramids in the moving and fixed images
    if( !m_MovingImagePyramid )
    {
        itkExceptionMacro(<<"Moving image pyramid is not present");
    }
    if( m_FixedMultiImagePyramid.size() != FImgTotal )
    {
        itkExceptionMacro(<<"The number of fixed images' pyramids must be equal to the number of fixed images.");
    }
    for( unsigned int f=0; f<FImgTotal; f++ )
    {
        if( !m_FixedMultiImagePyramid[f] )
        {
            itkExceptionMacro(<<"Fixed image pyramid " << f << " is not present" );
        }
    }

    // Are the fixed images' regions set? If not, use the buffered regions.
    if( this->m_FixedMultiImageRegion.size() == 0 )
    {
        for( unsigned int f=0; f<FImgTotal; f++ )
        {
            this->m_FixedMultiImageRegion.push_back( this->m_FixedMultiImage[f]->GetBufferedRegion() );
        }
    }
    else
    {
        if( this->m_FixedMultiImageRegion.size() != FImgTotal )
        {
            itkExceptionMacro( << "The number of regions must be equal to the "
                               << "number of fixed images" );
        }
    }

    // Setup the fixed image pyramid
    if( m_NumberOfLevelsSpecified )
    {
        for( unsigned int i=0; i<FImgTotal; i++ )
        {
            m_FixedMultiImagePyramid[i]->SetNumberOfLevels( m_NumberOfLevels );
        }
        m_MovingImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
    }

    if( m_ScheduleSpecified )
    {
        for( unsigned int i=0; i<FImgTotal; i++ )
        {
            m_FixedMultiImagePyramid[i]->SetNumberOfLevels( m_FixedImagePyramidSchedule.rows());
            m_FixedMultiImagePyramid[i]->SetSchedule( m_FixedImagePyramidSchedule );
        }

        m_MovingImagePyramid->SetNumberOfLevels( m_MovingImagePyramidSchedule.rows());
        m_MovingImagePyramid->SetSchedule( m_MovingImagePyramidSchedule );
    }

    for( unsigned int i=0; i<FImgTotal; i++ )
    {
        m_FixedMultiImagePyramid[i]->SetInput( this->m_FixedMultiImage[i] );
        m_FixedMultiImagePyramid[i]->UpdateLargestPossibleRegion();

        //el interpolador no necesita de piramide solo la region y esta es la misma
        //aplicada como interpolador actual de la imgen fija y no tiene nivel de piramide
        this->m_MultiInterpolator.push_back(this->GetMultiInterpolator()[i]);

    }

    // Setup the moving image pyramid
    m_MovingImagePyramid->SetInput( this->m_MovingImage );
    m_MovingImagePyramid->UpdateLargestPossibleRegion();

    typedef typename FixedImageRegionType::SizeType         SizeType;
    typedef typename FixedImageRegionType::IndexType        IndexType;

    ScheduleType movingschedule = m_MovingImagePyramid->GetSchedule();
    //itkDebugMacro ( << "MovingImage schedule: " << movingschedule );


    // Setup the pyramids for all the fixed images
    for( unsigned int f=0; f<FImgTotal; f++ )
    {
        ScheduleType schedule = this->m_FixedMultiImagePyramid[f]->GetSchedule();
        //itkDebugMacro ( << "FixedImage schedule: " << schedule );

        SizeType  inputSize  = this->m_FixedMultiImageRegion[f].GetSize();
        IndexType inputStart = this->m_FixedMultiImageRegion[f].GetIndex();

        const unsigned long NumberOfLevels =
                this->m_FixedMultiImagePyramid[f]->GetNumberOfLevels();

        this->m_FixedMultiImageRegionPyramid.push_back( FixedImageRegionPyramidType( NumberOfLevels ) );

        // Compute the FixedRegion corresponding to each level of the
        // pyramid. This uses the same algorithm of the ShrinkImageFilter
        // since the regions should be compatible.
        for ( unsigned int level=0; level < NumberOfLevels; level++ )
        {
            SizeType  size;
            IndexType start;
            for ( unsigned int dim = 0; dim < TFixedImage::ImageDimension; dim++)
            {
                const float scaleFactor = static_cast<float>( schedule[ level ][ dim ] );

                size[ dim ] = static_cast<typename SizeType::SizeValueType>(
                            vcl_floor(static_cast<float>( inputSize[ dim ] ) / scaleFactor ) );
                if( size[ dim ] < 1 )
                {
                    size[ dim ] = 1;
                }

                start[ dim ] = static_cast<typename IndexType::IndexValueType>(
                            vcl_ceil(static_cast<float>( inputStart[ dim ] ) / scaleFactor ) );
            }
            this->m_FixedMultiImageRegionPyramid[f][level].SetSize( size );
            this->m_FixedMultiImageRegionPyramid[f][level].SetIndex( start );
        }
    } // for( unsigned int f=0; f<FImgTotal; f++ )
}


/*
 * PrintSelf
 */
template < typename TFixedImage, typename TMovingImage >
void
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf( os, indent );

    for( unsigned int f=0; f<m_FixedMultiImagePyramid.size(); f++ )
    {
        os << indent << "FixedMultiImagePyramid[" << f <<"]: " ;
        os << m_FixedMultiImagePyramid[f].GetPointer() << std::endl;
    }
    os << indent << "MovingImagePyramid: ";
    os << m_MovingImagePyramid.GetPointer() << std::endl;

    os << indent << "NumberOfLevels: ";
    os << m_NumberOfLevels << std::endl;

    os << indent << "CurrentLevel: ";
    os << m_CurrentLevel << std::endl;

    os << indent << "InitialTransformParametersOfNextLevel: ";
    os << m_InitialTransformParametersOfNextLevel << std::endl;

    os << indent << "FixedImagePyramidSchedule : " << std::endl;
    os << m_FixedImagePyramidSchedule << std::endl;
    os << indent << "MovingImagePyramidSchedule : " << std::endl;
    os << m_MovingImagePyramidSchedule << std::endl;
}


/*
 * Generate Data
 */
template < typename TFixedImage, typename TMovingImage >
void
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::GenerateData()
{
    m_Stop = false;

    this->PreparePyramids();

    for ( m_CurrentLevel = 0; m_CurrentLevel < m_NumberOfLevels; m_CurrentLevel++ )
    {
        // Invoke an iteration event.
        // This allows a UI to reset any of the components between
        // resolution level.
        this->InvokeEvent( IterationEvent() );

        // Check if there has been a stop request
        if ( m_Stop )
        {
            break;
        }

        try
        {
            // Initialize the interconnects between components
            this->Initialize();
        }
        catch( ExceptionObject& err )
        {
            this->m_LastTransformParameters = ParametersType(1);
            this->m_LastTransformParameters.Fill( 0.0f );

            throw err;
        }

        // Do the optimization on each resolution level
        try
        {
            this->m_Optimizer->StartOptimization();
        }
        catch( ExceptionObject& err )
        {
            // An error has occurred in the optimization.
            // Update the parameters
            this->m_LastTransformParameters = this->m_Optimizer->GetCurrentPosition();

            // Pass exception to caller
            throw err;
        }

        // Get the results
        this->m_LastTransformParameters = this->m_Optimizer->GetCurrentPosition();
        this->m_Transform->SetParameters( this->m_LastTransformParameters );

        // Setup the initial parameters for next level
        if ( m_CurrentLevel < m_NumberOfLevels - 1 )
        {
            m_InitialTransformParametersOfNextLevel = this->m_LastTransformParameters;
        }
    } //for ( m_CurrentLevel = 0; m_CurrentLevel < m_NumberOfLevels; m_CurrentLevel++ )
}


/*
 * GetMTime
 */
template < typename TFixedImage, typename TMovingImage >
unsigned long
MultiResolutionMultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
::GetMTime() const
{
    unsigned long mtime = Superclass::GetMTime();
    unsigned long m;


    // Some of the following should be removed once ivars are put in the
    // input and output lists

    if ( this->m_Transform )
    {
        m = this->m_Transform->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    unsigned int InterpTotal = this->m_MultiInterpolator.size();
    for( unsigned int i=0; i<InterpTotal; i++ )
    {
        m = this->m_MultiInterpolator[i]->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    if ( this->m_MultiMetric )
    {
        m = this->m_MultiMetric->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    if ( this->m_Optimizer )
    {
        m = this->m_Optimizer->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    const unsigned int FImgTotal = this->GetNumberOfFixedImages();
    for( unsigned int i=0; i<FImgTotal; i++ )
    {
        m = this->m_FixedMultiImage[i]->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    if ( this->m_MovingImage )
    {
        m = this->m_MovingImage->GetMTime();
        mtime = (m > mtime ? m : mtime);
    }

    return mtime;

}

} // end namespace itk


#endif
