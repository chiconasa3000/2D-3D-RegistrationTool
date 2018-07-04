
#ifndef __itkNormalizedGradientCorrelationImageToImageMetric_h
#define __itkNormalizedGradientCorrelationImageToImageMetric_h

#include <itkImageToImageMetric.h>
#include <itkSobelOperator.h>
#include <itkNeighborhoodOperatorImageFilter.h>
#include <itkResampleImageFilter.h>


namespace itk
{


/** \class NormalizedGradientCorrelationImageToImageMetric
 * \brief Computes the similarity between two images using the Normalized Gradient Correlation metric.
 *
 * This class computes the similarity between two images by calculation of the
 * normalised correlation between their gradients. Sobel filters are used to
 * calculate the images' gradients along the fixed image's dimensions. Then,
 * individual calculations are made for the correlations in each dimension,
 * which are finally summed (not averaged).
 *
 * As the Normalized Gradient Correlation metric estimates the images
 * similarity using the images' gradient rather than their direct intensities,
 * it is less sensitive to low-frequency signals such as the ones introduced
 * by soft tissue.
 *
 * The class should work for 2 and 3 dimensional images. In case that the
 * fixed image has 3 dimensions but a single slice, the number of used
 * dimensions will be reduced accordingly.
 *
 * References: G. P. Penney et al. A comparison of similarity measures for use
 * in 2-D-3-D medical image registration. IEEE Transactions on Medical Imaging.
 * 17(4):586-95. August 1998.
 *
 * \warning Please note that the original implementation of the metric
 * averages the calculations made on each dimension, but this implementation
 * only sums them.
 *
 */
template < class TFixedImage, class TMovingImage >
class ITK_EXPORT NormalizedGradientCorrelationImageToImageMetric :
public ImageToImageMetric< TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef NormalizedGradientCorrelationImageToImageMetric Self;
  typedef ImageToImageMetric<TFixedImage, TMovingImage >  Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedGradientCorrelationImageToImageMetric, ImageToImageMetric);

/** Types transferred from the base class */
    typedef typename Superclass::RealType RealType;
// Work around a Visual Studio .NET bug
//#if defined(_MSC_VER) && (_MSC_VER == 1300)
//  typedef double RealType;
//#else
//  typedef typename Superclass::RealType                       RealType;
//#endif
  typedef typename Superclass::MovingImageType								MovingImageType;
  typedef typename Superclass::MovingImagePixelType						MovingImagePixelType;
  typedef typename Superclass::FixedImageType									FixedImageType;

  typedef typename FixedImageType::PixelType						      FixedImagePixelType;
  typedef ImageRegionConstIteratorWithIndex<FixedImageType>   FixedImageConstIteratorType;

  typedef typename Superclass::FixedImageConstPointer					FixedImageConstPointer;
  typedef typename Superclass::FixedImageRegionType						FixedImageRegionType;
  typedef typename FixedImageRegionType::SizeType                       SizeType;
  typedef typename Superclass::TransformType									TransformType;
  typedef typename Superclass::TransformPointer								TransformPointer;
  typedef typename Superclass::TransformParametersType				TransformParametersType;
  typedef typename Superclass::TransformJacobianType					TransformJacobianType;

  typedef  typename Superclass::InterpolatorType    InterpolatorType;
  typedef  typename Superclass::InterpolatorPointer InterpolatorPointer;

  typedef typename Superclass::MeasureType      MeasureType;
  typedef typename Superclass::DerivativeType   DerivativeType;

  typedef typename Superclass::InputPointType   InputPointType;

  typedef itk::ResampleImageFilter< MovingImageType,
                                    FixedImageType >  ResampleImageFilterType;

  itkStaticConstMacro(FixedImageDimension, unsigned int,
                                           TFixedImage::ImageDimension);

  typedef itk::Image<RealType,
                     itkGetStaticConstMacro(FixedImageDimension)> GradientImageType;

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
                            DerivativeType  & derivative ) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& Value, DerivativeType& derivative ) const;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject );

  /** Set/Get the value of Delta used for computing derivatives by finite
   * differences in the GetDerivative() method */
  itkSetMacro( DerivativeDelta, double );
  itkGetConstReferenceMacro( DerivativeDelta, double );

  /** Set the parameters defining the Transform. */
  void SetTransformParameters( const TransformParametersType & parameters ) const;

protected:
  NormalizedGradientCorrelationImageToImageMetric();
  virtual ~NormalizedGradientCorrelationImageToImageMetric() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  typedef NeighborhoodOperatorImageFilter<FixedImageType,
                                          GradientImageType> SobelFilterType;

private:
  NormalizedGradientCorrelationImageToImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&);                                  //purposely not implemented

  double  m_DerivativeDelta;
  unsigned int m_MaxDimension;

  /** The filter for transforming the moving images. */
  typename ResampleImageFilterType::Pointer m_ResampleImageFilter;

  /** The Sobel gradients of the fixed image */
  SobelOperator<RealType,
                itkGetStaticConstMacro(FixedImageDimension)>
                m_SobelOperators[FixedImageDimension];

  typename SobelFilterType::Pointer m_FixedSobelFilters[
    itkGetStaticConstMacro( FixedImageDimension ) ];
  typename SobelFilterType::Pointer m_MovingSobelFilters[
    itkGetStaticConstMacro( FixedImageDimension ) ];

  ZeroFluxNeumannBoundaryCondition<FixedImageType> m_FixedBoundaryCondition;
  ZeroFluxNeumannBoundaryCondition<FixedImageType> m_MovingBoundaryCondition;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNormalizedGradientCorrelationImageToImageMetric.txx"
#endif

#endif
