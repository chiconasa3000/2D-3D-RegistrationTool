#ifndef __itkPatternIntensityImageToImageMetric_txx
#define __itkPatternIntensityImageToImageMetric_txx


#include "itkPatternIntensityImageToImageMetric.h"
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkNumericTraits.h>


#define MaxMacro(a,b) ( a > b ? a : b )
#define MinMacro(a,b) ( a < b ? a : b )


namespace itk
{

/**
 * Constructor
 */
  template <class TFixedImage, class TMovingImage>
  PatternIntensityImageToImageMetric<TFixedImage,TMovingImage>
  ::PatternIntensityImageToImageMetric()
  {
    m_ResampleImageFilter = 0;
	  m_DifferenceImageFilter = 0;
    m_DerivativeDelta = 0.001;
	  m_Sigma = 10;
	  m_Radius = 3;
    m_MaxDimension = FixedImageDimension;
  }


/**
 * Initialize
 */
  template <class TFixedImage, class TMovingImage>
  void
  PatternIntensityImageToImageMetric<TFixedImage,TMovingImage>
  ::Initialize(void) throw ( ExceptionObject )
  {

    Superclass::Initialize();

    FixedImageRegionType fixedImageRegion =
      this->GetFixedImageRegion();
    typename FixedImageType::SizeType fixedImageRegionSize
      = fixedImageRegion.GetSize();
    typename FixedImageType::IndexType fixedImageRegionIndex
      = fixedImageRegion.GetIndex();

    typename FixedImageType::SizeType fixedImageSize =
      this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    typename FixedImageType::IndexType iterationStartIndex;
    typename FixedImageType::IndexType iterationEndIndex;
    typename FixedImageType::SizeType iterationSize;

    iterationStartIndex[0] = MaxMacro( fixedImageRegionIndex[0], m_Radius );
    iterationStartIndex[1] = MaxMacro( fixedImageRegionIndex[1], m_Radius );
    if( fixedImageRegionSize[2] > 1 ) // true for 3D images with more than one slice
      {
      iterationStartIndex[2] = MaxMacro( fixedImageRegionIndex[2], m_Radius );
      }
    else
      {
      iterationStartIndex[2] = 0;
      }

    iterationEndIndex[0] = MinMacro( fixedImageRegionIndex[0] + fixedImageRegionSize[0],
      fixedImageSize[0] - m_Radius );
    iterationEndIndex[1] = MinMacro( fixedImageRegionIndex[1] + fixedImageRegionSize[1],
      fixedImageSize[1] - m_Radius );
    if( fixedImageRegionSize[2] > 1 ) // true for 3D images with more than one slice
      {
      iterationEndIndex[2] = MinMacro( fixedImageRegionIndex[2] + fixedImageRegionSize[2],
        fixedImageSize[2] - m_Radius );
      }
    else
      {
      iterationEndIndex[2] = 1;
      }

    iterationSize[0] = iterationEndIndex[0] - iterationStartIndex[0];
    iterationSize[1] = iterationEndIndex[1] - iterationStartIndex[1];
    iterationSize[2] = iterationEndIndex[2] - iterationStartIndex[2];

    m_IterationRegion.SetIndex(iterationStartIndex);
    m_IterationRegion.SetSize(iterationSize);

    this->m_ResampleImageFilter = ResampleImageFilterType::New();
    this->m_ResampleImageFilter->SetInput( this->m_MovingImage );
    this->m_ResampleImageFilter->SetSize( this->m_FixedImage->GetLargestPossibleRegion().GetSize() );
    this->m_ResampleImageFilter->SetOutputSpacing( this->m_FixedImage->GetSpacing() );
    this->m_ResampleImageFilter->SetOutputOrigin( this->m_FixedImage->GetOrigin() );
    this->m_ResampleImageFilter->SetOutputDirection( this->m_FixedImage->GetDirection() );
    this->m_ResampleImageFilter->SetInterpolator( this->m_Interpolator );
    this->m_ResampleImageFilter->SetTransform( this->m_Transform );
    this->m_ResampleImageFilter->SetDefaultPixelValue( 0 );

    this->m_DifferenceImageFilter = DifferenceImageFilterType::New();
    this->m_DifferenceImageFilter->SetInput1( this->m_FixedImage );
    this->m_DifferenceImageFilter->SetInput2( this->m_ResampleImageFilter->GetOutput() );

    this->m_DifferenceImageFilter->GetOutput()->SetRequestedRegion( this->GetFixedImageRegion() );
  }


/**
 * PrintSelf
 */
  template <class TFixedImage, class TMovingImage>
  void
  PatternIntensityImageToImageMetric<TFixedImage,TMovingImage>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf( os, indent );
    os << indent << "DerivativeDelta: " << this->m_DerivativeDelta << std::endl;
    os << indent << "Radius ( in voxels ):" << this->m_Radius << std::endl;
    os << indent << "Sigma: " << this->m_Sigma << std::endl;
  }

/**
 * Get the value of the similarity measure
 */
  template <class TFixedImage, class TMovingImage>
  typename PatternIntensityImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
  PatternIntensityImageToImageMetric<TFixedImage,TMovingImage>
  ::GetValue( const TransformParametersType & parameters ) const
  {
    this->SetTransformParameters( parameters );

    try
      {
      this->m_DifferenceImageFilter->Update();
      }
    catch( ... )
      {
      itkExceptionMacro(<<
        "Calculation of difference image failed for PatternIntensity metric");
      }
    FixedImagePointer diffImage = this->m_DifferenceImageFilter->GetOutput();

    FixedImageConstIteratorType iter( diffImage, this->m_IterationRegion );

    MeasureType val = NumericTraits<MeasureType>::Zero;
    MeasureType diff = NumericTraits<MeasureType>::Zero;

    this->m_NumberOfPixelsCounted = 0;
    double SigmaSquared = this->m_Sigma*this->m_Sigma;
    for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
      {
      typename FixedImageType::IndexType fixedIndex = iter.GetIndex();

//Check if point is within the fixed image's mask
      InputPointType inputPoint;
      this->m_FixedImage->TransformIndexToPhysicalPoint( fixedIndex, inputPoint );

      if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        continue;
        }

      typename FixedImageType::IndexType kernelIndex = fixedIndex;

      int SliceRadius = ( this->m_IterationRegion.GetSize()[2] > 1 ?
        m_Radius : 0 );
      for( int k = -SliceRadius; k <= SliceRadius; k++ )
        {
        for(int j = -m_Radius; j<=m_Radius; j++)
          {
		      for(int i = -m_Radius; i<=m_Radius; i++)
            {
            kernelIndex[0] = fixedIndex[0] + i;
	          kernelIndex[1] = fixedIndex[1] + j;
            kernelIndex[2] = fixedIndex[2] + k;

            diff = diffImage->GetPixel( fixedIndex ) - diffImage->GetPixel( kernelIndex );
	          val += ( SigmaSquared )/( SigmaSquared + diff*diff );
	          }
          }
        }
      this->m_NumberOfPixelsCounted++;

      } //for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )

    if( this->m_NumberOfPixelsCounted == 0 )
      {
      itkExceptionMacro(<< "No voxels counted for metric calculation" );
      }

    //val /= this->m_NumberOfPixelsCounted;
    return val;
  }


/**
 * Get the Derivative Measure
 */
  template < class TFixedImage, class TMovingImage>
  void
  PatternIntensityImageToImageMetric<TFixedImage,TMovingImage>
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
  PatternIntensityImageToImageMetric<TFixedImage,TMovingImage>
  ::GetValueAndDerivative(const TransformParametersType & parameters,
                          MeasureType & Value, DerivativeType  & Derivative) const
  {
    Value      = this->GetValue( parameters );
    this->GetDerivative( parameters, Derivative );
  }

} // end namespace itk


#endif
