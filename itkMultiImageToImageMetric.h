#ifndef __itkMultiImageToImageMetric_h
#define __itkMultiImageToImageMetric_h

#include <itkArray.h>
#include <itkImageBase.h>
#include <itkInterpolateImageFunction.h>
#include <itkTransform.h>
#include <itkImageMaskSpatialObject.h>
#include <itkImageToImageMetric.h>
#include <itkSingleValuedCostFunction.h>
#include <itkTimeStamp.h>
#include <vector>

#include "itkMultiImageRegistrationMacro.h"

/** Macro used by the subclasses of itkMultiImageToImageMetric for creation of
 * individual metric objects. */
#define itkNewSingleMetricMacro(type) \
  virtual typename Superclass::MetricPointer NewSingleMetric (void) \
  { \
    return static_cast<type*> (type::New().GetPointer()); \
  }


namespace itk
{

/** \class MultiImageToImageMetric
 * \brief Computes similarity between a moving image and multiple fixed images.
 *
 * This class computes the similarity between a moving image and multiple
 * fixed images as the sum of the individual similarity measures between the
 * moving image and each of the fixed ones.
 *
 * This class is abstract and its subclasses are expected to provide concrete
 * implementation for different similarity metrics, such as mean squares or
 * pattern intensity. In many cases, only a header (.h) file is required to
 * implement a metric, which should set the MetricType typedef accordingly
 * (see itkMeanSquaresMultiImageToImageMetric.h and
 * itkNormalizedCorrelationMultiImageToImageMetric for examples)
 *
 * The fixed images are connected to this class as a std::vector of images.
 * Also, the fixed images' regions and interpolators must come in
 * vectors, which must have the same number of elements. Like ITK's
 * ImageToImageMetric, this class expects a single transform to be plugged in.
 *
 */
  template < class TFixedImage, class TMovingImage >
  class ITK_EXPORT MultiImageToImageMetric : public SingleValuedCostFunction
  {
  public:
/** Standard class typedefs */
    typedef MultiImageToImageMetric         Self;
    typedef SingleValuedCostFunction        Superclass;
    typedef SmartPointer<Self>              Pointer;
    typedef SmartPointer<const Self>        ConstPointer;

/** Type used for representing point components  */
    typedef typename Superclass::ParametersValueType CoordinateRepresentationType;

/** Run-time type information (and related methods). */
    itkTypeMacro(MultiImageToImageMetric, SingleValuedCostFunction);

/**  Type of the moving Image. */
    typedef TMovingImage                               MovingImageType;
    typedef typename TMovingImage::PixelType           MovingImagePixelType;
    typedef typename MovingImageType::ConstPointer     MovingImageConstPointer;

/**  Type of the fixed Images. */
    typedef TFixedImage                             FixedImageType;
    typedef typename FixedImageType::Pointer        FixedImagePointer;
    typedef typename FixedImageType::ConstPointer   FixedImageConstPointerType;
    typedef std::vector<FixedImageConstPointerType> FixedMultiImageType;

    typedef typename FixedImageType::RegionType     FixedImageRegionType;
    typedef std::vector<FixedImageRegionType>       FixedMultiImageRegionType;

/** Constants for the image dimensions */
    itkStaticConstMacro(MovingImageDimension, unsigned int, MovingImageType::ImageDimension);
    itkStaticConstMacro(FixedImageDimension, unsigned int, FixedImageType::ImageDimension);

/**  Type of the Transform Base class */
    typedef Transform<CoordinateRepresentationType,
      itkGetStaticConstMacro(MovingImageDimension),
      itkGetStaticConstMacro(FixedImageDimension)>  TransformType;
    typedef typename TransformType::Pointer         TransformPointer;
    typedef typename TransformType::InputPointType  InputPointType;
    typedef typename TransformType::OutputPointType OutputPointType;
    typedef typename TransformType::ParametersType  TransformParametersType;
    typedef typename TransformType::JacobianType    TransformJacobianType;

/**  Type of the Interpolator Base class */
    typedef InterpolateImageFunction<MovingImageType,
      CoordinateRepresentationType>             InterpolatorType;
    typedef typename InterpolatorType::Pointer   InterpolatorPointer;
    typedef std::vector<InterpolatorPointer>    MultiInterpolatorType;

///** Gaussian filter to compute the gradient of the Moving Image */
//    typedef typename NumericTraits<MovingImagePixelType>::RealType
//      RealType;
//    typedef CovariantVector<RealType,
//      itkGetStaticConstMacro(MovingImageDimension)> GradientPixelType;
//    typedef Image<GradientPixelType,
//      itkGetStaticConstMacro(MovingImageDimension)> GradientImageType;
//    typedef SmartPointer<GradientImageType> GradientImagePointer;
//    typedef GradientRecursiveGaussianImageFilter< MovingImageType,
//      GradientImageType > GradientImageFilterType;
//    typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;

/**  Type for the mask of the fixed image. Only pixels that are "inside"
     this mask will be considered for the computation of the metric */
    typedef ImageMaskSpatialObject< itkGetStaticConstMacro(FixedImageDimension) >
      FixedImageMaskType;
    typedef typename FixedImageMaskType::Pointer       FixedImageMaskPointer;
    typedef typename FixedImageMaskType::ConstPointer  FixedImageMaskConstPointer;
    typedef std::vector<FixedImageMaskPointer>    FixedMultiImageMaskType;

///**  Type for the mask of the moving image. Only pixels that are "inside"
//     this mask will be considered for the computation of the metric */
//    typedef SpatialObject< itkGetStaticConstMacro(MovingImageDimension) >
//      MovingImageMaskType;
//    typedef typename MovingImageMaskType::Pointer      MovingImageMaskPointer;
//    typedef typename MovingImageMaskType::ConstPointer MovingImageMaskConstPointer;

/**  Type of the measure. */
    typedef typename Superclass::MeasureType MeasureType;

/**  Type of the derivative. */
    typedef typename Superclass::DerivativeType DerivativeType;

/**  Type of the parameters. */
    typedef typename Superclass::ParametersType ParametersType;

/** Type of the single-input metrics. */
    typedef ImageToImageMetric<TFixedImage, TMovingImage>  MetricType;
    typedef typename MetricType::Pointer                   MetricPointer;
    typedef std::vector<MetricPointer>                     MultiMetricType;

/** Set or get the Fixed Images.  */
    itkSetConstStdVectorMacro(FixedMultiImage,FixedMultiImageType);
    itkGetConstStdVectorMacro(FixedMultiImage,FixedMultiImageType);

/** Set or get the Moving Image.  */
    itkSetConstObjectMacro( MovingImage, MovingImageType );
    itkGetConstObjectMacro( MovingImage, MovingImageType );

/** Set or get the Transform. */
    itkSetObjectMacro( Transform, TransformType );
    itkGetConstObjectMacro( Transform, TransformType );

/** Add, set or get the Interpolators. */
    itkAddObjectToStdVectorMacro(Interpolator,InterpolatorType,MultiInterpolator);
    itkSetConstStdVectorMacro(MultiInterpolator,MultiInterpolatorType);
    itkGetConstStdVectorMacro(MultiInterpolator,MultiInterpolatorType);

/** Set or get the region over which the metric will be computed */
    itkSetConstStdVectorMacro(FixedMultiImageRegion,FixedMultiImageRegionType);
    itkGetConstStdVectorMacro(FixedMultiImageRegion,FixedMultiImageRegionType);

///** Set/Get the moving image mask. */
//  itkSetObjectMacro( MovingImageMask, MovingImageMaskType );
//  itkSetConstObjectMacro( MovingImageMask, MovingImageMaskType );
//  itkGetConstObjectMacro( MovingImageMask, MovingImageMaskType );
//

/** Set or Get the fixed images' masks. */
    itkSetConstStdVectorMacro(FixedMultiImageMask,FixedMultiImageMaskType);
    itkGetConstStdVectorMacro(FixedMultiImageMask,FixedMultiImageMaskType);

/** Set ot Get the step size used for derivative calculation */
    itkSetMacro( DerivativeDelta, double);
    itkGetConstReferenceMacro( DerivativeDelta, double );

///** Set or Get gradient computation. */
//    itkSetMacro( ComputeGradient, bool);
//    itkGetConstReferenceMacro( ComputeGradient, bool);
//    itkBooleanMacro(ComputeGradient);
//
///** Computes the gradient image and assigns it to m_GradientImage */
//    virtual void ComputeGradient();
//
///** Get Gradient Image. */
//    itkGetConstObjectMacro( GradientImage, GradientImageType );

/** Set the parameters defining the Transform. */
    void SetTransformParameters( const ParametersType & parameters ) const;

/** Return the number of parameters required by the Transform */
    unsigned int GetNumberOfParameters(void) const
      { return m_Transform->GetNumberOfParameters(); }

/** Initialize the Metric by making sure that all the components
  *  are present and plugged together correctly     */
    virtual void Initialize(void) throw ( ExceptionObject );

    virtual void DoNumberRevision(void) const throw ( ExceptionObject );
    virtual void DoConnectionRevision(void) const throw ( ExceptionObject );
    virtual void DoFullRevision(void) const throw ( ExceptionObject )
      { DoNumberRevision(); DoConnectionRevision();}

/** Get the derivatives of the match measure. */
    virtual void GetDerivative( const TransformParametersType & parameters,
      DerivativeType & derivative ) const;

/**  Get the value for single valued optimizers. */
    virtual MeasureType GetValue( const TransformParametersType & parameters ) const;

  protected:
    MultiImageToImageMetric();
    virtual ~MultiImageToImageMetric();
    void PrintSelf(std::ostream& os, Indent indent) const;

/** Create a new instance of an individual metric. */
    virtual MetricPointer NewSingleMetric()=0;

    FixedMultiImageType         m_FixedMultiImage;
    MovingImageConstPointer     m_MovingImage;

    mutable TransformPointer    m_Transform;
    MultiInterpolatorType       m_MultiInterpolator;

    //bool                        m_ComputeGradient;
    //GradientImagePointer        m_GradientImage;

    FixedMultiImageMaskType     m_FixedMultiImageMask;
    //MovingImageMaskConstPointer m_MovingImageMask;

    MultiMetricType             m_MultiMetric;

    double                      m_DerivativeDelta;

  private:
    MultiImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&);          //purposely not implemented

    FixedMultiImageRegionType m_FixedMultiImageRegion;

    TimeStamp m_InitializationTime;
  };

} //namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiImageToImageMetric.txx"
#endif

#endif
