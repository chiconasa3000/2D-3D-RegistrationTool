
#ifndef __itkGradientDifferenceSingleImageToImageMetric_h
#define __itkGradientDifferenceSingleImageToImageMetric_h

#include <itkImageToImageMetric.h>

#include <itkSobelOperator.h>
#include <itkNeighborhoodOperatorImageFilter.h>
#include <itkResampleImageFilter.h>

namespace itk
{

/** \class GradientDifferenceSingleImageToImageMetric
 * \brief Computes the similarity between two images using the Gradient Difference metric.
 *
 * The Gradient Difference metric estimates the similarity between two images
 * based on the the image obtained by subtraction between the gradients of the
 * fixed and moving images. Gradient difference voxels considered for
 * calculation are summed using  a reciprocal form, which should make the
 * metric strong against outliers present in the images, such as surgical
 * instruments. Besides, taking the gradients rather than the direct
 * intensities reduces the effects of low frequency signals, such as the ones
 * introduced by soft tissue.
 *
 * This metric can be used on 3D and 2D images. In case of 3D fixed images
 * with a single slice, only two dimensions will be considered for
 * the metric calculation.
 *
 * This metric does not differ much from the original
 * itk::GradientDifferenceImageToImageMetric class. The main difference is
 * that the original ITK class precalculates the moving image's gradient and
 * this implementation does not. As this class was designed to be used with
 * the itk::GradientDifferenceMultiImageToImageMetric class, keeping the
 * gradient precalculation would be an innecessary waste of resources, as each
 * of the individual metrics would calculate the gradient and store it,
 * keeping multiple copies of the same image in memory.
 *
 * References: G. P. Penney et al. A comparison of similarity measures for use
 * in 2-D-3-D medical image registration. IEEE Transactions on Medical Imaging.
 * 17(4):586-95. August 1998.
 *
 * \warning Gradient Difference can be defined with a factor 's' that multiplies
 * the moving image's gradient before subtraction. This factor cannot be set in
 * the current implementation and is always equal to 1.
 *
 */
template < class TFixedImage, class TMovingImage >
class ITK_EXPORT GradientDifferenceSingleImageToImageMetric :
public ImageToImageMetric< TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef GradientDifferenceSingleImageToImageMetric           Self;
  typedef ImageToImageMetric<TFixedImage, TMovingImage > Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GradientDifferenceSingleImageToImageMetric, ImageToImageMetric);


/** Types transferred from the base class */
// Work around a Visual Studio .NET bug
#if defined(_MSC_VER) && (_MSC_VER == 1300)
  typedef double RealType;
#else
  typedef typename Superclass::RealType                 RealType;
#endif
  typedef typename Superclass::TransformType            TransformType;
  typedef typename Superclass::TransformPointer         TransformPointer;
  typedef typename Superclass::TransformParametersType  TransformParametersType;
  typedef typename Superclass::TransformJacobianType    TransformJacobianType;

  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
  typedef typename Superclass::MovingImageConstPointer  MovingImageConstPointer;
  typedef typename Superclass::InputPointType           InputPointType;

  typedef typename FixedImageType::RegionType           FixedImageRegionType;
  typedef typename FixedImageRegionType::IndexType      IndexType;
  typedef typename FixedImageRegionType::SizeType       SizeType;

  typedef typename TFixedImage::PixelType               FixedImagePixelType;
  typedef typename TMovingImage::PixelType              MovedImagePixelType;

  itkStaticConstMacro(FixedImageDimension, unsigned int, TFixedImage::ImageDimension);
  /** Types for transforming the moving image */

  typedef itk::ResampleImageFilter< MovingImageType,
                                    FixedImageType >  ResampleImageFilterType;

  typedef RealType                               GradientPixelType;
  typedef itk::Image<GradientPixelType,
            itkGetStaticConstMacro(FixedImageDimension)>  GradientImageType;

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

protected:
  GradientDifferenceSingleImageToImageMetric();
  virtual ~GradientDifferenceSingleImageToImageMetric() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Compute the variance and range of the moving image gradients. */
  void ComputeVariance( void ) const;

  typedef NeighborhoodOperatorImageFilter<
    FixedImageType, GradientImageType > SobelFilterType;

private:
  GradientDifferenceSingleImageToImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The variance of the moving image gradients. */
  mutable GradientPixelType m_Variance[FixedImageDimension];

  /** The filter for transforming the moving image. */
  typename ResampleImageFilterType::Pointer m_ResampleImageFilter;

  SobelOperator< GradientPixelType,
                 itkGetStaticConstMacro(FixedImageDimension) >
               m_SobelOperators[FixedImageDimension];

  typename SobelFilterType::Pointer m_FixedSobelFilter[
    itkGetStaticConstMacro( FixedImageDimension )];
  typename SobelFilterType::Pointer m_MovingSobelFilter[
    itkGetStaticConstMacro( FixedImageDimension )];

  ZeroFluxNeumannBoundaryCondition<FixedImageType >
    m_MovingBoundaryCondition;
  ZeroFluxNeumannBoundaryCondition<FixedImageType >
    m_FixedBoundaryCondition;

  double        m_DerivativeDelta;
  unsigned int  m_MaxDimension;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGradientDifferenceSingleImageToImageMetric.txx"
#endif

#endif
