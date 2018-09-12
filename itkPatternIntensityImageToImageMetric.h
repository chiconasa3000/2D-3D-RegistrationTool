#ifndef __itkPatternIntensityImageToImageMetric_h
#define __itkPatternIntensityImageToImageMetric_h


#include <itkImageToImageMetric.h>
#include <itkResampleImageFilter.h>
#include <itkSubtractImageFilter.h>


namespace itk
{
/** \class PatternIntensityImageToImageMetric
 * \brief Computes Pattern intensity similarity metric between two images to be registered
 *
 * The Pattern Intensity metric calculates the similarity between two images
 * based on the intensities of the image obtained by subtraction of the fixed
 * and moving images. A moving kernel is passed over this image detecting
 * 'structures', defined as the areas with large intensity variation. Voxel
 * values are summed using a reciprocal form, which gives low values to areas
 * with structures and large ones to areas without.
 *
 * This metric, differing from the Gradiend Difference and Normalized Gradient
 * Correlation metrics, relies on the images' direct intensity values rather
 * than their gradients, but its moving kernel makes it insensitive to
 * low-frequency variations introduced by soft tissue. Its reciprocal form
 * makes it robust against outliers such as surgical instruments.
 *
 * Pattern Intensity has two parameters with are the Radius and Sigma. Radius
 * defines the size of the moving kernel in voxels, which is a cube of
 * 2*Radius + 1 voxels per side. The Sigma parameter defines the metric's
 * sensitivity to structure. As a rule of thumb, it should have a value
 * similar to the fixed and moving images' standard deviation (which implies
 * that fixed and moving images should have comparable intensity levels).
 *
 * The used implementation was taken from Penney, G., et al. (1998), "A
 * comparison of similarity metrics for use un 2D-3D medical image
 * registration", IEEE Transactions in Medical Imaging, 17(4):586-595.
 *
 * \warning This implementation of the Pattern Intensity metric only works
 * for images with three dimensions (or images with three dimensions and a
 * single slice)
 *
 */
template < class TFixedImage, class TMovingImage >
class ITK_EXPORT PatternIntensityImageToImageMetric :
public ImageToImageMetric< TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef PatternIntensityImageToImageMetric           Self;
  typedef ImageToImageMetric<TFixedImage, TMovingImage > Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PatternIntensityImageToImageMetric, ImageToImageMetric);


  /** Typedefs from the superclass. */
  typedef typename
    Superclass::CoordinateRepresentationType              CoordinateRepresentationType;
  typedef typename Superclass::MovingImageType            MovingImageType;
  typedef typename Superclass::MovingImagePixelType       MovingImagePixelType;
  typedef typename Superclass::MovingImageType::Pointer	  MovingImagePointer;
  typedef typename Superclass::MovingImageConstPointer    MovingImageConstPointer;

  typedef typename Superclass::FixedImageType             FixedImageType;
  typedef typename FixedImageType::Pointer                FixedImagePointer;
  typedef typename Superclass::FixedImageConstPointer     FixedImageConstPointer;
  typedef typename Superclass::FixedImageRegionType       FixedImageRegionType;

  typedef typename Superclass::TransformType              TransformType;
  typedef typename Superclass::TransformPointer           TransformPointer;
  typedef typename Superclass::InputPointType             InputPointType;
  typedef typename Superclass::OutputPointType            OutputPointType;
  typedef typename Superclass::TransformParametersType    TransformParametersType;
  typedef typename Superclass::TransformJacobianType      TransformJacobianType;

  typedef typename Superclass::InterpolatorType           InterpolatorType;
  typedef typename Superclass::InterpolatorPointer        InterpolatorPointer;

  typedef typename Superclass::RealType                   RealType;
  typedef typename Superclass::GradientPixelType          GradientPixelType;
  typedef typename Superclass::GradientImageType          GradientImageType;
  typedef typename Superclass::GradientImagePointer       GradientImagePointer;
  typedef typename Superclass::GradientImageFilterType    GradientImageFilterType;
  typedef typename Superclass::GradientImageFilterPointer GradientImageFilterPointer;
  typedef typename Superclass::FixedImageMaskType         FixedImageMaskType;
  typedef typename Superclass::FixedImageMaskPointer      FixedImageMaskPointer;
  typedef typename Superclass::MovingImageMaskType        MovingImageMaskType;
  typedef typename Superclass::MovingImageMaskPointer     MovingImageMaskPointer;
  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::ParametersType             ParametersType;

  typedef typename Superclass::FixedImageType::PixelType    FixedImagePixelType;
  typedef typename Superclass::MovingImageType::RegionType  MovingImageRegionType;

  typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedImageConstIteratorType;

  typedef itk::ResampleImageFilter< MovingImageType,
                                    FixedImageType >  ResampleImageFilterType;

  typedef itk::SubtractImageFilter< FixedImageType,
                                    FixedImageType,
                                    FixedImageType > DifferenceImageFilterType;

  /** The fixed image dimension. */
  itkStaticConstMacro( FixedImageDimension, unsigned int,
    FixedImageType::ImageDimension );

  /** The moving image dimension. */
  itkStaticConstMacro( MovingImageDimension, unsigned int,
    MovingImageType::ImageDimension );

  /** Get the value for single valued optimizers. */
  virtual MeasureType GetValue( const TransformParametersType & parameters ) const;

  /** Get the derivatives of the match measure. */
  virtual void GetDerivative( const TransformParametersType & parameters,
    DerivativeType & derivative ) const;

  /** Get value and derivatives for multiple valued optimizers. */
  virtual void GetValueAndDerivative( const TransformParametersType & parameters,
    MeasureType& Value, DerivativeType& Derivative ) const;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly.*/
  virtual void Initialize(void) throw ( ExceptionObject );

/** Set or Get the sigma constant  */
  itkSetMacro( Sigma , double );
  itkGetConstReferenceMacro( Sigma, double );

/** Set or Get the radius  */
  itkSetMacro( Radius , int );
  itkGetConstReferenceMacro( Radius, int );

protected:
  PatternIntensityImageToImageMetric();
  virtual ~PatternIntensityImageToImageMetric() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  PatternIntensityImageToImageMetric(const Self&);  //purposely not implemented
  void operator=(const Self&);                      //purposely not implemented

	typename ResampleImageFilterType::Pointer   m_ResampleImageFilter;
	typename DifferenceImageFilterType::Pointer m_DifferenceImageFilter;

	double  m_DerivativeDelta;
	int     m_Radius;
	double  m_Sigma;
  unsigned int m_MaxDimension;

	FixedImageRegionType	m_IterationRegion;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPatternIntensityImageToImageMetric.txx"
#endif

#endif
