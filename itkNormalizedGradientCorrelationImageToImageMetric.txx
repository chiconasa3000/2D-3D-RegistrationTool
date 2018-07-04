
#ifndef __itkNormalizedGradientCorrelationImageToImageMetric_txx
#define __itkNormalizedGradientCorrelationImageToImageMetric_txx

#include "itkNormalizedGradientCorrelationImageToImageMetric.h"
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkNumericTraits.h>


namespace itk
{

/**
 * Constructor
 */
  template <class TFixedImage, class TMovingImage>
  NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>
  ::NormalizedGradientCorrelationImageToImageMetric()
  {
    m_DerivativeDelta = 0.001;
    m_MaxDimension = FixedImageDimension;
  }


/**
 * Initialize
 */
  template <class TFixedImage, class TMovingImage>
  void
  NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>
  ::Initialize(void) throw ( ExceptionObject )
  {
// Initialise the base class
    Superclass::Initialize();

    SizeType size = this->GetFixedImageRegion().GetSize();
    for( unsigned int dim=0; dim < FixedImageDimension; dim++ )
      {
      if( size[dim] < 2 )
        {
        m_MaxDimension = dim;
        break;
        }
      }

    for (unsigned int dim=0; dim<m_MaxDimension; dim++)
      {
      m_SobelOperators[dim].SetRadius( 1 );
      m_SobelOperators[dim].SetDirection( dim );
      m_SobelOperators[dim].CreateDirectional();

      m_FixedSobelFilters[dim] = SobelFilterType::New();
      m_FixedSobelFilters[dim]->OverrideBoundaryCondition(
        &m_FixedBoundaryCondition );
      m_FixedSobelFilters[dim]->SetOperator( m_SobelOperators[dim] );
      m_FixedSobelFilters[dim]->SetInput( this->GetFixedImage() );
      m_FixedSobelFilters[dim]->GetOutput()->SetRequestedRegion(
        this->GetFixedImageRegion() );
      }

    m_ResampleImageFilter = ResampleImageFilterType::New();

    m_ResampleImageFilter->SetTransform( this->m_Transform );
    m_ResampleImageFilter->SetInterpolator( this->m_Interpolator );
    m_ResampleImageFilter->SetInput( this->m_MovingImage );
    m_ResampleImageFilter->SetDefaultPixelValue(
      itk::NumericTraits<FixedImagePixelType>::Zero );
    m_ResampleImageFilter->UseReferenceImageOn();
    m_ResampleImageFilter->SetReferenceImage( this->m_FixedImage );
    m_ResampleImageFilter->GetOutput()->SetRequestedRegion(
      this->GetFixedImageRegion() );
    m_ResampleImageFilter->Update();

    for (unsigned int dim=0; dim < m_MaxDimension; dim++)
      {
      m_MovingSobelFilters[dim] = SobelFilterType::New();
      m_MovingSobelFilters[dim]->OverrideBoundaryCondition(
        &m_MovingBoundaryCondition );
      m_MovingSobelFilters[dim]->SetOperator( m_SobelOperators[dim] );
      m_MovingSobelFilters[dim]->SetInput(m_ResampleImageFilter->GetOutput());
      m_MovingSobelFilters[dim]->GetOutput()->SetRequestedRegion(
        this->GetFixedImageRegion() );
      }
  }


/**
 * PrintSelf
 */
  template <class TFixedImage, class TMovingImage>
  void
  NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf( os, indent );
    os << indent << "DerivativeDelta: " << this->m_DerivativeDelta << std::endl;
  }


/**
 * Get the value of the similarity measure
 */
  template <class TFixedImage, class TMovingImage>
  typename NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
  NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>
  ::GetValue( const TransformParametersType & parameters ) const
  {
    this->m_NumberOfPixelsCounted = 0;
    this->SetTransformParameters( parameters );

    for (unsigned int dim=0; dim<m_MaxDimension; dim++)
      {
      this->m_FixedSobelFilters[dim]->Update();
      this->m_MovingSobelFilters[dim]->Update();
      }

    MeasureType val = NumericTraits< MeasureType >::Zero;

/*
  cc: cross corrrelation
  fac: fixed auto correlation, this is, auto correlation of the fixed image
  mac: moving auto correlation, this is, moving image auto correlation
*/
    MeasureType cc[FixedImageDimension];
    MeasureType fac[FixedImageDimension];
    MeasureType mac[FixedImageDimension];

    for( unsigned int dim=0; dim < m_MaxDimension; dim++ )
      {
      cc[dim] = NumericTraits< MeasureType >::Zero;
      fac[dim] = NumericTraits< MeasureType >::Zero;
      mac[dim] = NumericTraits< MeasureType >::Zero;
      }

    RealType movingGradient[FixedImageDimension];
    RealType fixedGradient[FixedImageDimension];

    FixedImageConstIteratorType iter( this->m_FixedImage,
      this->GetFixedImageRegion() );
    for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      typename FixedImageType::IndexType fixedIndex = iter.GetIndex();

//Check if point is inside the fixed image mask
      InputPointType inputPoint;
      this->GetFixedImage()->TransformIndexToPhysicalPoint( fixedIndex, inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        continue;
        }

      for( unsigned int dim=0; dim<m_MaxDimension; dim++ )
        {
        fixedGradient[dim] = m_FixedSobelFilters[dim]->GetOutput()->GetPixel(
          fixedIndex );
        movingGradient[dim] = m_MovingSobelFilters[dim]->GetOutput()->GetPixel(
          fixedIndex );
        cc[dim] +=  movingGradient[dim] * fixedGradient[dim];
        fac[dim] += fixedGradient[dim] * fixedGradient[dim];
        mac[dim] += movingGradient[dim] * movingGradient[dim];
        }

      this->m_NumberOfPixelsCounted++;
      }

    if( this->m_NumberOfPixelsCounted == 0 )
      {
      itkExceptionMacro(<< "No voxels counted for metric calculation");
      }

    for( unsigned int dim=0; dim < m_MaxDimension; dim++ )
      {
      if( fac[dim] == NumericTraits< MeasureType >::Zero ||
          mac[dim] == NumericTraits< MeasureType >::Zero )
        {
        itkExceptionMacro(<< "Auto correlation(s) equal to zero");
        }
      }

    for( unsigned int dim=0; dim < m_MaxDimension; dim++ )
      {
      val += cc[dim] / vcl_sqrt( fac[dim] * mac[dim] ) / m_MaxDimension;
      }
    return val;
  }


/**
 * Get the Derivative Measure
 */
  template < class TFixedImage, class TMovingImage>
  void
  NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>
  ::GetDerivative( const TransformParametersType & parameters,
                   DerivativeType & derivative           ) const
  {
    TransformParametersType testPoint;
    testPoint = parameters;

    const unsigned int numberOfParameters = this->GetNumberOfParameters();
    derivative = DerivativeType( numberOfParameters );

    for( unsigned int i=0; i<numberOfParameters; i++)
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
 * Get both the match Measure and theDerivative Measure
 */
  template <class TFixedImage, class TMovingImage>
  void
  NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>
  ::GetValueAndDerivative(const TransformParametersType & parameters,
                          MeasureType & Value, DerivativeType  & Derivative) const
  {
    Value      = this->GetValue( parameters );
    this->GetDerivative( parameters, Derivative );
  }


/**
 * Set the parameters that define a unique transform
 */
  template <class TFixedImage, class TMovingImage>
  void
  NormalizedGradientCorrelationImageToImageMetric<TFixedImage,TMovingImage>
  ::SetTransformParameters( const TransformParametersType & parameters ) const
  {
    if( !this->m_Transform )
      {
      itkExceptionMacro(<<"Transform has not been assigned");
      }
    this->m_Transform->SetParameters( parameters );
  }

} // end namespace itk

#endif
