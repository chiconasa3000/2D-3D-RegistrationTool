#ifndef __itkMultiImageToImageRegistrationMethod_h
#define __itkMultiImageToImageRegistrationMethod_h

#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkSingleValuedNonLinearOptimizer.h>
#include <itkDataObjectDecorator.h>

#include "itkMultiImageToImageMetric.h"

/** \class MultiImageToImageRegistrationMethod
 * \brief Base class for registration methods with multiple fixed images.
 *
 * The MultiImageToImageRegistrationMethod class coordinates all objects
 * requiered in the registration process, which are the images, transform,
 * metric, optimizer and interpolators. As this class is a subclass of
 * itk::ProcessObject, a call to the Update method is all required to
 * execute the registration after all components have been plugged in.
 *
 * The class is templated over the classes of the fixed and moving images.
 * The moving image, transform and optimizer are set using the SetMovingImage,
 * SetTransform and SetOptimizer methods respectively. As in the original
 * ITK class, the optimizer must be a single valued cost function.
 *
 * The class expects a number of fixed images for registration, which are
 * stored internally on a std::vector object. It is possible to give the
 * array of fixed images to the class using the SetFixedMultiImage method
 * or add them one by one using the AddFixedImage method. The second way
 * is recommended, as it leads to simpler code.
 *
 * Like the fixed images, the MultiImageToImageRegistrationMethod class
 * requieres multiple interpolators, regions and masks. Again, it is
 * recommended to add these objects into the registration method using the
 * corresponding Add method, which are AddInterpolator, AddFixedImageRegion
 * and AddFixedImageMask. In general terms, the number of interpolators,
 * regions and masks must be equal to the number of fixed images. However, if
 * the fixed images' regions are not given, their buffered regions will be
 * used by default. The fixed images' masks can also be left unset, which will
 * make them be ignored during registration. Also, it is possible to define
 * masks for only some of the fixed images. To do this, call SetFixedImageMask
 * giving NULL as argument for the images that should not be masked.
 *
 * The metric must be a subclass of MultiImageToImageMetric, which returns the
 * sum of a set of individual metrics calculated between the moving image and
 * each of the fixed ones, using the corresponding interpolator and transform.
 *
 * \warning It is not possible to define masks for the moving image in the
 * current implementation.
 *
 */
namespace itk
{

template <typename TFixedImage, typename TMovingImage>
class ITK_EXPORT MultiImageToImageRegistrationMethod : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MultiImageToImageRegistrationMethod Self;
  typedef ProcessObject                       Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiImageToImageRegistrationMethod, ProcessObject);

  /**  Type of the Fixed image. */
  typedef          TFixedImage                    FixedImageType;
  typedef typename FixedImageType::Pointer        FixedImagePointer;
  typedef typename FixedImageType::ConstPointer   FixedImageConstPointer;
  typedef std::vector<FixedImageConstPointer>     FixedMultiImageType;

  /**  Type of the Moving image. */
  typedef          TMovingImage                     MovingImageType;
  typedef typename MovingImageType::ConstPointer    MovingImageConstPointer;

  /**  Type of the metric. */
  typedef MultiImageToImageMetric<FixedImageType,MovingImageType> MultiMetricType;
  typedef typename MultiMetricType::Pointer                       MultiMetricPointer;

  /** Type of the fixed images regions */
  typedef typename MultiMetricType::FixedImageRegionType          FixedImageRegionType;
  typedef typename MultiMetricType::FixedMultiImageRegionType     FixedMultiImageRegionType;

  /** Type of the fixed images' masks */
  typedef typename MultiMetricType::FixedImageMaskType          FixedImageMaskType;
  typedef typename MultiMetricType::FixedImageMaskConstPointer  FixedImageMaskConstPointer;
  typedef typename MultiMetricType::FixedMultiImageMaskType     FixedMultiImageMaskType;

  /**  Type of the Transform . */
  typedef  typename MultiMetricType::TransformType  TransformType;
  typedef  typename TransformType::Pointer          TransformPointer;

  /** Type for the output: Using Decorator pattern for enabling
   *  the Transform to be passed in the data pipeline */
  typedef  DataObjectDecorator< TransformType >      TransformOutputType;
  typedef typename TransformOutputType::Pointer      TransformOutputPointer;
  typedef typename TransformOutputType::ConstPointer TransformOutputConstPointer;

  /**  Type of the Interpolator. */
  typedef  typename MultiMetricType::InterpolatorType       InterpolatorType;
  typedef  typename InterpolatorType::Pointer               InterpolatorPointer;
  typedef  typename MultiMetricType::MultiInterpolatorType  MultiInterpolatorType;

  /**  Type of the optimizer. */
  typedef SingleValuedNonLinearOptimizer  OptimizerType;
  typedef typename OptimizerType::Pointer OptimizerPointer;

  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef  typename MultiMetricType::TransformParametersType  ParametersType;

  /** Smart Pointer type to a DataObject. */
  typedef typename DataObject::Pointer DataObjectPointer;

  /** Set/Get the Fixed multi-image. */
  virtual void SetFixedMultiImage( const FixedMultiImageType );
  itkGetConstStdVectorMacro( FixedMultiImage, FixedMultiImageType );
  virtual void AddFixedImage( const FixedImageType* );

  /** Set or Get the Moving image. */
  virtual void SetMovingImage( const MovingImageType * movingImage );
  itkGetConstObjectMacro( MovingImage, MovingImageType );

  /** Set or Get the Optimizer. */
  itkSetObjectMacro( Optimizer,  OptimizerType );
  itkGetObjectMacro( Optimizer,  OptimizerType );

  /** Set or Get the MultiMetric */
  itkSetObjectMacro( MultiMetric, MultiMetricType );
  itkGetObjectMacro( MultiMetric, MultiMetricType );

  /** Set or Get the Transfrom. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetObjectMacro( Transform, TransformType );

  /** Set or Get the MultiInterpolator. */
  itkSetConstStdVectorMacro( MultiInterpolator, MultiInterpolatorType );
  itkGetConstStdVectorMacro( MultiInterpolator, MultiInterpolatorType );
  itkAddObjectToStdVectorMacro( Interpolator, InterpolatorType, MultiInterpolator );

  /** Set ot Get the fixed images' regions */
  itkSetConstStdVectorMacro( FixedMultiImageRegion, FixedMultiImageRegionType );
  itkGetConstStdVectorMacro( FixedMultiImageRegion, FixedMultiImageRegionType );
  itkAddToStdVectorMacro( FixedImageRegion, FixedImageRegionType, FixedMultiImageRegion );

  /** Add, Set or Get the fixed images' masks */
  itkAddObjectToStdVectorMacro( FixedImageMask, FixedImageMaskType, FixedMultiImageMask );
  itkGetConstStdVectorMacro( FixedMultiImageMask, FixedMultiImageMaskType );
  itkSetConstStdVectorMacro( FixedMultiImageMask, FixedMultiImageMaskType );

  /** Set or Get the initial transformation parameters. */
  virtual void SetInitialTransformParameters( const ParametersType & param );
  itkGetConstReferenceMacro( InitialTransformParameters, ParametersType );

  /** Get the last transformation parameters visited by
   * the optimizer. */
  itkGetConstReferenceMacro( LastTransformParameters, ParametersType );

  /** Get the number of fixed images present in the registration */
  const unsigned int GetNumberOfFixedImages() const;

  ///** Set the region of the fixed image to be considered as region of
  // interest during the registration. This region will be passed to
  // the ImageMetric in order to restrict the metric computation to
  // consider only this region.
  // \warning The same region can also be set directly into the metric.
  // please avoid to set the region in both places since this can lead
  // to inconsistent configurations.  */
  //void SetFixedImageRegion( const  FixedImageRegionType & region );
  ///** Get the region of the fixed image to be considered as region of
  // interest during the registration. This region will be passed to
  // the ImageMetric in order to restrict the metric computation to
  // consider only this region.  */
  //itkGetConstReferenceMacro( FixedImageRegion, FixedImageRegionType );
  ///** True if a region has been defined for the fixed image to which
  // the ImageMetric will limit its computation */
  //itkGetConstMacro( FixedImageRegionDefined, bool );
  ///** Turn on/off the use of a fixed image region to which
  // the ImageMetric will limit its computation.
  // \warning The region must have been previously defined using the
  // SetFixedImageRegion member function */
  //itkSetMacro( FixedImageRegionDefined, bool );

  /** Initialize by setting the interconnects between the components. */
  void Initialize() throw (ExceptionObject);

  /** Returns the transform resulting from the registration process  */
  virtual const TransformOutputType * GetOutput() const;

  /** Make a DataObject of the correct type to be used as the specified
   * output. */
  DataObjectPointer MakeOutput(unsigned int idx);

  /** Method to return the latest modified time of this object or
   * any of its cached ivars */
  unsigned long GetMTime() const;

protected:
  MultiImageToImageRegistrationMethod();
  virtual ~MultiImageToImageRegistrationMethod() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  void  GenerateData ();

  /** Provides derived classes with the ability to set this private var */
  itkSetMacro( LastTransformParameters, ParametersType );

  /** Method that initiates the optimization process. This method should not be
   * called directly by the users. Instead, this method is intended to be
   * invoked internally by the Update() method. */
  void StartOptimization(void);

  MovingImageConstPointer   m_MovingImage;
  FixedMultiImageType       m_FixedMultiImage;

  TransformPointer          m_Transform;
  MultiMetricPointer        m_MultiMetric;
  OptimizerType::Pointer    m_Optimizer;
  MultiInterpolatorType     m_MultiInterpolator;

  ParametersType            m_InitialTransformParameters;
  ParametersType            m_LastTransformParameters;

  FixedMultiImageRegionType m_FixedMultiImageRegion;
  FixedMultiImageMaskType   m_FixedMultiImageMask;


private:
  MultiImageToImageRegistrationMethod(const Self&); //purposely not implemented
  void operator=(const Self&);                      //purposely not implemented
};


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiImageToImageRegistrationMethod.txx"
#endif

#endif
