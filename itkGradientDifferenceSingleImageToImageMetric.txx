#ifndef __itkGradientDifferenceSingleImageToImageMetric_txx
#define __itkGradientDifferenceSingleImageToImageMetric_txx

#include "itkGradientDifferenceSingleImageToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNumericTraits.h"

namespace itk
{

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage>
GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>
::GradientDifferenceSingleImageToImageMetric()
{
  m_ResampleImageFilter = 0;

  m_MaxDimension = FixedImageDimension;

  for (unsigned dim=0; dim < m_MaxDimension; dim++)
    {
    m_Variance[dim] = itk::NumericTraits<GradientPixelType>::Zero;
    }

  this->m_DerivativeDelta = 0.001;
}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage>
void
GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
{
// Initialise the base class
  Superclass::Initialize();

  SizeType size = this->GetFixedImageRegion().GetSize();
  for (unsigned int dim=0; dim<FixedImageDimension; dim++)
    {
    if( size[dim] < 2 )
      {
      m_MaxDimension = dim;
      break;
      }
    }

  // Create the filter to resample the moving image

  m_ResampleImageFilter = ResampleImageFilterType::New();

  m_ResampleImageFilter->SetTransform(    this->m_Transform );
  m_ResampleImageFilter->SetInterpolator( this->m_Interpolator );
  m_ResampleImageFilter->SetInput( this->m_MovingImage );

  m_ResampleImageFilter->SetDefaultPixelValue( 0 );

  m_ResampleImageFilter->SetSize( this->m_FixedImage->GetLargestPossibleRegion().GetSize() );
  m_ResampleImageFilter->SetOutputOrigin( this->m_FixedImage->GetOrigin() );
  m_ResampleImageFilter->SetOutputSpacing( this->m_FixedImage->GetSpacing() );
  m_ResampleImageFilter->SetOutputDirection( this->m_FixedImage->GetDirection() );

  m_ResampleImageFilter->GetOutput()->SetRequestedRegion( this->GetFixedImageRegion() );

  // Compute the gradient of the fixed image
  for (unsigned int dim=0; dim<m_MaxDimension; dim++)
    {
    m_SobelOperators[dim].SetRadius( 1 );
    m_SobelOperators[dim].SetDirection( dim );
    m_SobelOperators[dim].CreateDirectional();

    m_FixedSobelFilter[dim] = SobelFilterType::New();

    m_FixedSobelFilter[dim]->OverrideBoundaryCondition( &m_FixedBoundaryCondition );
    m_FixedSobelFilter[dim]->SetOperator( m_SobelOperators[dim] );

    m_FixedSobelFilter[dim]->SetInput( this->m_FixedImage );

    m_FixedSobelFilter[dim]->GetOutput()->SetRequestedRegion( this->GetFixedImageRegion() );
    m_FixedSobelFilter[dim]->Update();
    }

  this->ComputeVariance();

  // Compute the gradient of the transformed moving image
  for (unsigned int dim=0; dim<m_MaxDimension; dim++)
    {
    m_MovingSobelFilter[dim] = SobelFilterType::New();

    m_MovingSobelFilter[dim]->OverrideBoundaryCondition( &m_MovingBoundaryCondition );
    m_MovingSobelFilter[dim]->SetOperator( m_SobelOperators[dim] );

    m_MovingSobelFilter[dim]->SetInput( m_ResampleImageFilter->GetOutput() );

    m_MovingSobelFilter[dim]->GetOutput()->SetRequestedRegion( this->GetFixedImageRegion() );
    m_MovingSobelFilter[dim]->Update();
    }

}


/**
 * PrintSelf
 */
template <class TFixedImage, class TMovingImage>
void
GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "DerivativeDelta: " << this->m_DerivativeDelta << std::endl;
}


/**
 * Compute the gradient variances in each dimension.
 */
template <class TFixedImage, class TMovingImage>
void
GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>
::ComputeVariance( void ) const
{
  unsigned long nPixels;
  GradientPixelType mean[FixedImageDimension];

  for (unsigned int dim=0; dim < m_MaxDimension; dim++)
    {

    typedef itk::ImageRegionConstIteratorWithIndex<
      GradientImageType > IteratorType;

    IteratorType iterate( m_FixedSobelFilter[dim]->GetOutput(),
                          this->GetFixedImageRegion() );

// Calculate the gradients' mean first, and their range.
    nPixels =  0;

    mean[dim] = itk::NumericTraits<GradientPixelType>::Zero;

    for ( iterate.GoToBegin(); !iterate.IsAtEnd(); ++iterate )
      {
      InputPointType inputPoint;
      this->m_FixedImage->TransformIndexToPhysicalPoint( iterate.GetIndex(), inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        continue;
        }

      mean[dim] += iterate.Get();
      nPixels++;
      }

    if (nPixels > 1)
      {
      mean[dim] /= nPixels;
      }
      else
      {
      itkExceptionMacro(<<"Insufficient pixels for mean and variance calculations");
      }

// Now calculate the variance
    m_Variance[dim] = itk::NumericTraits<GradientPixelType>::Zero;
    GradientPixelType gradient;

    for ( iterate.GoToBegin(); !iterate.IsAtEnd(); ++iterate )
      {
      InputPointType inputPoint;
      this->m_FixedImage->TransformIndexToPhysicalPoint( iterate.GetIndex(), inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        continue;
        }

      gradient = iterate.Get();
      gradient -= mean[dim];

      m_Variance[dim] += gradient*gradient;
      }

    m_Variance[dim] /= ( nPixels - 1 );

    } //for (unsigned int dim=0; dim < m_MaxDimension; dim++)
}


/**
 * Get the value of the similarity measure
 */
template <class TFixedImage, class TMovingImage>
typename GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>
::GetValue( const TransformParametersType & parameters ) const
{
  this->SetTransformParameters( parameters );

  m_ResampleImageFilter->Update();

  for (unsigned int dim=0; dim < m_MaxDimension; dim++)
    {
    m_MovingSobelFilter[dim]->Update();
    }

  typedef  itk::ImageRegionConstIteratorWithIndex< FixedImageType >
    IteratorType;

  MeasureType val = NumericTraits< MeasureType >::Zero;
  this->m_NumberOfPixelsCounted = 0;

  IteratorType iter( this->m_FixedImage, this->GetFixedImageRegion() );
  for ( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    IndexType index = iter.GetIndex();

//Check if the fixed image's point is inside the mask
    InputPointType inputPoint;
    this->m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

    if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
      {
      continue;
      }


    for( unsigned int dim=0; dim<m_MaxDimension; dim++ )
      {
      if( m_Variance[dim] == itk::NumericTraits<GradientPixelType>::Zero )
        {
        continue;
        }

      GradientPixelType diff =
        m_FixedSobelFilter[dim]->GetOutput()->GetPixel( index ) -
        m_MovingSobelFilter[dim]->GetOutput()->GetPixel( index );

        val += m_Variance[dim] / ( m_Variance[dim] + diff * diff );
      }

    this->m_NumberOfPixelsCounted++;
    }

  if( this->m_NumberOfPixelsCounted == 0 )
    {
    itkExceptionMacro(<< "No pixels counted for metric calculation" );
    }

  return val;
}


/**
 * Get the Derivative Measure
 */
template < class TFixedImage, class TMovingImage>
void
GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>
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
    testPoint[i] += 2 * this->m_DerivativeDelta;
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
GradientDifferenceSingleImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters,
                        MeasureType & Value, DerivativeType  & Derivative) const
{
  Value      = this->GetValue( parameters );
  this->GetDerivative( parameters, Derivative );
}

} // end namespace itk


#endif
