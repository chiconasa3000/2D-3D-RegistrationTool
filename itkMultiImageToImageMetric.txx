#ifndef __itkMultiImageToImageMetric_txx
#define __itkMultiImageToImageMetric_txx

#include "itkMultiImageToImageMetric.h"

namespace itk {

/**
 * Constructor
 */
  template <class TFixedImage, class TMovingImage>
  MultiImageToImageMetric<TFixedImage,TMovingImage>
    ::MultiImageToImageMetric()
  {
    m_MovingImage       = 0; // has to be provided by the user.
    m_Transform         = 0; // has to be provided by the user.
    //m_GradientImage   = 0; // will receive the output of the filter;
    //m_ComputeGradient = false;
    m_DerivativeDelta = 0.001;
  }

/**
 * Destructor
 */
  template <class TFixedImage, class TMovingImage>
  MultiImageToImageMetric<TFixedImage,TMovingImage>::
  ~MultiImageToImageMetric()
  {
    const unsigned int MetricTotal = m_MultiMetric.size();
    for( unsigned int metricNum=0; metricNum<MetricTotal; metricNum++)
    {
      m_MultiMetric.back() = NULL;
      m_MultiMetric.pop_back();
    }
  }


/**
 * DoConnectionRevision
 */
  template <class TFixedImage, class TMovingImage>
  void
  MultiImageToImageMetric<TFixedImage,TMovingImage>
  ::DoConnectionRevision( ) const throw ( ExceptionObject )
  {
    const unsigned int FImgTotal = m_FixedMultiImage.size();
    for( unsigned int fImg=0; fImg < FImgTotal; fImg++ )
    {
      if( m_MultiMetric[fImg]->GetFixedImage() != m_FixedMultiImage[fImg] )
        {
        itkExceptionMacro( << "Incorrect connection between metric " << fImg << " and its fixed image" );
        }

      if( m_MultiMetric[fImg]->GetFixedImageRegion() != m_FixedMultiImageRegion[fImg] )
        {
        itkExceptionMacro( << "Incorrect connection between metric " << fImg << " and its region" );
        }

      if( m_FixedMultiImageMask.size() > 0 && m_FixedMultiImageMask[fImg] &&
        m_MultiMetric[fImg]->GetFixedImageMask() != m_FixedMultiImageMask[fImg] )
        {
        itkExceptionMacro( << "Incorrect connection between metrics " << fImg << " and its mask" );
        }

      if( m_MultiMetric[fImg]->GetMovingImage() != m_MovingImage )
        {
        itkExceptionMacro( << "Incorrect connection between metric " << fImg << " and the moving image" );
        }

      //if( m_MultiMetric[fImg]->GetMovingImageMask() != m_MovingImageMask )
      //  {
      //    itkExceptionMacro( << "Incorrect connection between metric " << fImg << " and the moving image's mask" );
      //  }

      if( m_MultiMetric[fImg]->GetInterpolator() != m_MultiInterpolator[fImg] )
        {
        itkExceptionMacro( << "Incorrect connection between metric " << fImg << " and its interpolator" );
        }

      if( m_MultiMetric[fImg]->GetTransform() != m_Transform )
        {
        itkExceptionMacro( << "Incorrect connection between metric " << fImg << " and the transform" );
        }
      }
  }

/**
 * DoNumberRevision
 */
  template <class TFixedImage, class TMovingImage>
  void
  MultiImageToImageMetric<TFixedImage,TMovingImage>
  ::DoNumberRevision( ) const throw ( ExceptionObject )
  {
//Does the moving image exist?
    if( !m_MovingImage )
      {
      itkExceptionMacro( << "Moving Image is not present" );
      }

//Does the transform exist?
    if( !m_Transform )
      {
      itkExceptionMacro( << "Transform is not present" );
      }

// Are there any fixed images?
    const unsigned int FImgTotal = m_FixedMultiImage.size();
    if( FImgTotal == 0 )
      {
      itkExceptionMacro( << "Fixed Images vector is empty" );
      }

//Does the interpolator vector have the same number of elements as the number
//of fixed images?
    const unsigned int InterpTotal = m_MultiInterpolator.size();
    if( InterpTotal != FImgTotal )
      {
      itkExceptionMacro( << "Interpolators vector does not have the same number "
        "of elements as the total fixed images" );
      }

//Now do the usual checks on the regions...
    const unsigned int RegionsTotal = m_FixedMultiImageRegion.size();
    if( RegionsTotal != FImgTotal )
      {
      itkExceptionMacro( <<"Fixed Images Regions vector does not have the "
        "same number of elements as the total fixed images" );
      }

//Number of metrics == Number of fixed images ?
    const unsigned int MetricTotal = m_MultiMetric.size();
    if( MetricTotal != FImgTotal )
      {
      itkExceptionMacro(<<"Metrics vector does not have the same number "
        "of elements as the total fixed images");
      }

//Do the fixed images' mask array have the correct number of elements?
//Remember can be empty or have the same number of elements as fixed images.
//Any other value is not allowed.
    const unsigned int FixedMaskTotal = m_FixedMultiImageMask.size();
    if( FixedMaskTotal > 0 && FixedMaskTotal != FImgTotal )
      {
      itkExceptionMacro(<<"Fixed images' mask array has an incorrect number "
        "of elements. Allowed values are 0 or the number  of fixed images");
      }

// Check arrays for any null elements
    for( unsigned int fImgNum=0; fImgNum < FImgTotal; fImgNum++ )
      {
      FixedImageConstPointerType fImg = m_FixedMultiImage[ fImgNum ];

//Do all fixed images exist?
      if( !fImg )
        {
        itkExceptionMacro( << "Fixed Images vector contains null element(s)" );
        }

//Check interpolators...
      InterpolatorPointer interp = m_MultiInterpolator[ fImgNum ];
      if( !interp )
        {
        itkExceptionMacro( << "Interpolators vector contains null element(s)" );
        }

      FixedImageRegionType fRegion = m_FixedMultiImageRegion[ fImgNum ];
      unsigned long pixelTotal = fRegion.GetNumberOfPixels();
      if( pixelTotal == 0 )
      {
        itkExceptionMacro(<<"Fixed image region is empty");
      }
      if ( !fRegion.Crop( fImg->GetBufferedRegion() ) )
      {
        itkExceptionMacro(
          <<"FixedImageRegion does not overlap the fixed image buffered region" );
      }

// Check metrics
      MetricPointer metric = m_MultiMetric[ fImgNum ];
      if( !metric )
        {
        itkExceptionMacro( << "Metrics vector contains null element(s)" );
        }

// No searches for null elements are made on the fixed images' masks array.
// It is valid that some elements are defined while others are not.

    } // for( int fImgNum=0; fImgNum < FImgTotal; fImgNum++ )

    if( m_DerivativeDelta <= 0.0 )
      {
      itkExceptionMacro( << "Derivative delta must be greater than zero" );
      }
  }


/**
 * PrintSelf
 */
  template <class TFixedImage, class TMovingImage>
  void
  MultiImageToImageMetric<TFixedImage,TMovingImage>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf( os, indent );
    os << indent << "Moving Image: " << m_MovingImage.GetPointer()
      << std::endl;
    for( unsigned int i=0; i<m_FixedMultiImage.size(); i++ )
      {
      os << indent << "FixedMultiImage[" << i <<"]: " << m_FixedMultiImage[i]
         << std::endl;
      }

    os << indent << "Transform:     " << m_Transform.GetPointer() << std::endl;
    for( unsigned int i=0; i<m_MultiInterpolator.size(); i++ )
      {
        os << indent << "MultiInterpolator[" << i << "]: " <<
          m_MultiInterpolator[i] << std::endl;
      }

    for( unsigned int i=0; i<m_MultiMetric.size(); i++ )
      {
        os << indent << "MultiMetric[" << i <<"]: " << m_MultiMetric[i] <<
          std::endl;
      }

    //os << indent << "ComputeGradient: "
    //   << static_cast<typename NumericTraits<bool>::PrintType>(m_ComputeGradient)
    //   << std::endl;
    //os << indent << "Gradient Image: " << m_GradientImage.GetPointer()
    //   << std::endl;

    for( unsigned int i=0; i<m_MultiMetric.size(); i++ )
      {
      os << indent << "FixedMultiImageRegion[" << i << "]: " <<
        m_FixedMultiImageRegion[i] << std::endl;
      }

    for( unsigned int i=0; i<m_FixedMultiImageMask.size(); i++ )
      {
      os << indent << "FixedMultiImageMask[" << i << "]: ";
      if( m_FixedMultiImageMask[i] )
        {
        os << m_FixedMultiImageMask[i].GetPointer() << std::endl;
        }
        else
        {
        os << "0" << std::endl;
        }
      }

    //os << indent << "Moving Image Mask: " << m_MovingImageMask.GetPointer()
    //   << std::endl;

    os << indent << "Derivative delta: " << m_DerivativeDelta << std::endl;
  }

/**
 * Set the parameters that define a unique transform
 */
  template <class TFixedImage, class TMovingImage>
  void
  MultiImageToImageMetric<TFixedImage,TMovingImage>
  ::SetTransformParameters( const ParametersType & parameters ) const
  {
    if( !m_Transform )
      {
      itkExceptionMacro(<<"Transform has not been assigned");
      }
    m_Transform->SetParameters( parameters );
  }

/**
 * Initialize
 */
  template <class TFixedImage, class TMovingImage>
  void
  MultiImageToImageMetric<TFixedImage,TMovingImage>
  ::Initialize(void) throw ( ExceptionObject )
  {
    const unsigned int FImgTotal = m_FixedMultiImage.size();
    if( FImgTotal == 0 )
      {
      itkExceptionMacro( << "Fixed Images vector is empty" );
      }

// Allocate the individual metrics
    const unsigned int OldMetricTotal = m_MultiMetric.size();
    for( unsigned int i=0; i<OldMetricTotal; i++ )
    {
      m_MultiMetric.back() = NULL;
      m_MultiMetric.pop_back(); // calls 'delete' on last element
    }

    for( unsigned int fImgNum=0; fImgNum<FImgTotal; fImgNum++ )
      {
      m_MultiMetric.push_back( NewSingleMetric() );
      }

// Now check the number of objects by calling DoNumberRevision()
    DoNumberRevision();

// If the images are provided by different sources, update them.
    if( m_MovingImage->GetSource() )
      {
      m_MovingImage->GetSource()->Update();
      }

    for( unsigned int fImgNum=0; fImgNum < FImgTotal; fImgNum++ )
      {
      FixedImageConstPointerType fImg = m_FixedMultiImage[ fImgNum ];

      if( fImg->GetSource() )
        {
        fImg->GetSource()->Update();
        }

      InterpolatorPointer interp = m_MultiInterpolator[fImgNum];
      interp->SetInputImage( m_MovingImage );
      }

////Calculate the Moving Image gradient, if necessary
//    if ( m_ComputeGradient )
//      {
//      this->ComputeGradient();
//      }

    for( unsigned int fImg=0; fImg < FImgTotal; fImg++ )
      {
      m_MultiMetric[fImg]->SetFixedImage( m_FixedMultiImage[fImg] );
      m_MultiMetric[fImg]->SetFixedImageRegion( m_FixedMultiImageRegion[fImg] );
      if( m_FixedMultiImageMask.size() == FImgTotal )
        {
        m_MultiMetric[fImg]->SetFixedImageMask( m_FixedMultiImageMask[fImg] );
        }
      m_MultiMetric[fImg]->SetMovingImage( m_MovingImage );
      //m_MultiMetric[fImg]->SetMovingImageMask( m_MovingImageMask );
      m_MultiMetric[fImg]->SetInterpolator( m_MultiInterpolator[fImg] );
      m_MultiMetric[fImg]->SetTransform( m_Transform );

//ComputeGradient is set to false on each of the individual metrics. We don't
//want to calculate the same gradient multiple times.
      m_MultiMetric[fImg]->SetComputeGradient( false );

      m_MultiMetric[fImg]->Initialize();
      }

// Update the object's initialization time
    m_InitializationTime.Modified();

// If there are any observers on the metric, call them to give the
// user code a chance to set parameters on the metric
    this->InvokeEvent( InitializeEvent() );
  }


///*
// * Compute the gradient image and assign it to m_GradientImage.
// */
//  template <class TFixedImage, class TMovingImage>
//  void
//  MultiImageToImageMetric<TFixedImage,TMovingImage>
//  ::ComputeGradient()
//  {
//    GradientImageFilterPointer gradientFilter
//      = GradientImageFilterType::New();
//
//    gradientFilter->SetInput( m_MovingImage );
//
//    const typename MovingImageType::SpacingType&
//      spacing = m_MovingImage->GetSpacing();
//    double maximumSpacing=0.0;
//    for(unsigned int i=0; i<MovingImageDimension; i++)
//      {
//      if( spacing[i] > maximumSpacing )
//        {
//        maximumSpacing = spacing[i];
//        }
//      }
//    gradientFilter->SetSigma( maximumSpacing );
//    gradientFilter->SetNormalizeAcrossScale( true );
//
//  #ifdef ITK_USE_ORIENTED_IMAGE_DIRECTION
//    gradientFilter->SetUseImageDirection( true );
//  #endif
//
//    gradientFilter->Update();
//
//    m_GradientImage = gradientFilter->GetOutput();
//  }

/**
 * Get the derivatives of the match measure.
 */
  template <class TFixedImage, class TMovingImage>
  void
  MultiImageToImageMetric<TFixedImage,TMovingImage>
  ::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
  {

  itkDebugMacro( "GetDerivative( " << parameters << ", " << derivative << " ) " );

  if( this->GetMTime() >= m_InitializationTime.GetMTime() )
    {
    this->DoFullRevision();
    }

  TransformParametersType testPoint;
  testPoint = parameters;

  const unsigned int ParTotal = this->GetNumberOfParameters();
  derivative = DerivativeType( ParTotal );

  for( unsigned int i=0; i<ParTotal; i++)
    {
    testPoint[i] -= this->m_DerivativeDelta;
    const MeasureType valuep0 = this->GetValue( testPoint );
    testPoint[i] += 2* this->m_DerivativeDelta;
    const MeasureType valuep1 = this->GetValue( testPoint );
    derivative[i] = (valuep1 - valuep0 ) / ( 2 * this->m_DerivativeDelta );

    testPoint[i] = parameters[i];
    }
  }

/**
 * Get the value for single valued optimizers.
 */
  template <class TFixedImage, class TMovingImage>
  typename MultiImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
  MultiImageToImageMetric<TFixedImage,TMovingImage>
  ::GetValue( const TransformParametersType & parameters ) const
  {
    itkDebugMacro( "GetValue( " << parameters << " ) " );

    if( this->GetMTime() >= m_InitializationTime.GetMTime() )
      {
      this->DoFullRevision();
      }

    this->SetTransformParameters( parameters );

    MeasureType measure = NumericTraits< MeasureType >::Zero;

    unsigned int FImgTotal = m_FixedMultiImage.size();
    for( unsigned int fImgNum=0; fImgNum<FImgTotal; fImgNum++ )
    {
      measure += m_MultiMetric[fImgNum]->GetValue( parameters );
    }
    return measure;
  }

} //namespace itk

#endif
