#ifndef __itkMultiResolutionMultiImageToImageRegistrationMethod_h
#define __itkMultiResolutionMultiImageToImageRegistrationMethod_h

#include <itkProcessObject.h>
#include <itkSingleValuedNonLinearOptimizer.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkNumericTraits.h>
#include <itkDataObjectDecorator.h>

#include "itkMultiImageRegistrationMacro.h"
#include "itkMultiImageToImageMetric.h"
#include "itkMultiImageToImageRegistrationMethod.h"


/** \class MultiResolutionMultiImageToImageRegistrationMethod
 * \brief Multi-resolution registration using multiple fixed images.
 *
 * The MultiResolutionMultiImageToImageRegistrationMethod class coordinates
 * all required objects for registration with multiple resolution levels.
 * Registration is repeated for each resolution level, giving the result of
 * each level as the starting point for the following one. A careful
 * selection of the number of levels and downsampling factors can considerably
 * increase the registration speed and robustness, avoiding convergence into
 * local minima.
 *
 * This class is a subclass of MultiImageToImageRegistrationMethod, so the
 * methods for definition of the fixed and moving images, interpolators,
 * transform, optimizer, regions and masks are present. In addition to these
 * components, multi resolution registration requires filters to generate the
 * downsampled versions of the moving and fixed images, known as pyramids. The
 * moving image's filter is set with the SetMovingImagePyramid method and the
 * fixed images' ones can be set with one call to SetFixedMultiImagePyramid or
 * added one by one using the AddFixedImagePyramid. Again, the second method
 * is recommended for its simplicity. In all cases, the number of pyramid
 * filters must be equal to the number of fixed images.
 *
 * The class offers two ways to define the number of resolution levels, which
 * are the SetNumberOfLevels of the SetSchedules methods. When using
 * SetNumberOfLevels only the total number of resolution levels must be given,
 * and the class will generate a default set of downsampling factors. Using
 * SetSchedules is more complex, as two matrices must be given -one for moving
 * and one for the fixed images- which define both the number of levels (rows)
 * and the downsampling factors applied to each dimension (columns). In the
 * case of 2D-3D registration, the use of SetSchedules is always recommended,
 * as the downsampling factor in the fixed images' slice direction must be
 * always equal to one (i.e. the fixed images' schedule must have its third
 * column filled with ones).
 *
 *
 */
namespace itk
{
template <typename TFixedImage, typename TMovingImage>
class ITK_EXPORT MultiResolutionMultiImageToImageRegistrationMethod
  : public MultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionMultiImageToImageRegistrationMethod             Self;
  typedef MultiImageToImageRegistrationMethod<TFixedImage,TMovingImage>  Superclass;
  typedef SmartPointer<Self>                                             Pointer;
  typedef SmartPointer<const Self>                                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionMultiImageToImageRegistrationMethod,
    MultiImageToImageRegistrationMethod);

  /**  Type of the Fixed images. */
  typedef typename Superclass::FixedImageType         FixedImageType;
  typedef typename Superclass::FixedImagePointer      FixedImagePointer;
  typedef typename Superclass::FixedImageConstPointer FixedImageConstPointer;
  typedef typename Superclass::FixedMultiImageType    FixedMultiImageType;

  /** Type of the fixed images' regions */
  typedef typename Superclass::FixedImageRegionType       FixedImageRegionType;
  typedef typename Superclass::FixedMultiImageRegionType  FixedMultiImageRegionType;
  typedef std::vector<FixedImageRegionType>               FixedImageRegionPyramidType;
  typedef std::vector<FixedImageRegionPyramidType>        FixedMultiImageRegionPyramidType;

  typedef typename Superclass::FixedImageMaskType         FixedImageMaskType;
  typedef typename Superclass::FixedImageMaskConstPointer FixedImageMaskConstPointer;
  typedef typename Superclass::FixedMultiImageMaskType    FixedMultiImageMaskType;

  /**  Type of the Moving image. */
  typedef typename Superclass::MovingImageType         MovingImageType;
  typedef typename Superclass::MovingImageConstPointer MovingImageConstPointer;

  /**  Type of the multi metric. */
  typedef typename Superclass::MultiMetricType    MultiMetricType;
  typedef typename Superclass::MultiMetricPointer MultiMetricPointer;

  /**  Type of the Transform */
  typedef typename Superclass::TransformType    TransformType;
  typedef typename Superclass::TransformPointer TransformPointer;

  /** Type for the output: Using Decorator pattern for enabling
   *  the Transform to be passed in the data pipeline */
  typedef typename Superclass::TransformOutputType          TransformOutputType;
  typedef typename Superclass::TransformOutputPointer       TransformOutputPointer;
  typedef typename Superclass::TransformOutputConstPointer  TransformOutputConstPointer;

  /**  Type of the Interpolator. */
  typedef typename Superclass::InterpolatorType       InterpolatorType;
  typedef typename Superclass::InterpolatorPointer    InterpolatorPointer;
  typedef typename Superclass::MultiInterpolatorType  MultiInterpolatorType;

  /**  Type of the optimizer. */
  typedef typename Superclass::OptimizerType    OptimizerType;
  typedef typename Superclass::OptimizerPointer OptimizerPointer;

  /** Type of the Fixed image multiresolution pyramid. */
  typedef MultiResolutionPyramidImageFilter<
    FixedImageType, FixedImageType >              FixedImagePyramidType;
  typedef typename FixedImagePyramidType::Pointer FixedImagePyramidPointer;
  typedef std::vector<FixedImagePyramidPointer>   FixedMultiImagePyramidType;

  /** Type of pyramid schedule type */
  typedef typename FixedImagePyramidType::ScheduleType ScheduleType;

  /** Type of the moving image multiresolution pyramid. */
  typedef MultiResolutionPyramidImageFilter<
    MovingImageType, MovingImageType >             MovingImagePyramidType;
  typedef typename MovingImagePyramidType::Pointer MovingImagePyramidPointer;

  /** Type of the Transformation parameters This is the same type used to
   *  represent the search space of the optimization algorithm */
  typedef typename Superclass::ParametersType ParametersType;

  /** Smart Pointer type to a DataObject. */
  typedef typename Superclass::DataObjectPointer DataObjectPointer;

  /** Method to stop the registration. */
  void StopRegistration();

  /** Set/Get the Fixed image pyramid. */
  itkSetConstStdVectorMacro( FixedMultiImagePyramid, FixedMultiImagePyramidType );
  itkGetConstStdVectorMacro( FixedMultiImagePyramid, FixedMultiImagePyramidType );
  itkAddObjectToStdVectorMacro( FixedImagePyramid, FixedImagePyramidType, FixedMultiImagePyramid );

  /** Set/Get the Moving image pyramid. */
  itkSetObjectMacro( MovingImagePyramid, MovingImagePyramidType );
  itkGetObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Set or get the schedules . */
  void SetSchedules( const ScheduleType & fixedSchedule,
                     const ScheduleType & movingSchedule );
  itkGetConstMacro( FixedImagePyramidSchedule, ScheduleType );
  itkGetConstMacro( MovingImagePyramidSchedule, ScheduleType );

  /** Set or get the number of multi-resolution levels. */
  void SetNumberOfLevels( unsigned long numberOfLevels );
  itkGetConstMacro( NumberOfLevels, unsigned long );

  /** Get the current resolution level being processed. */
  itkGetConstMacro( CurrentLevel, unsigned long );

  /** Set/Get the initial transformation parameters of the next resolution
   level to be processed. The default is the last set of parameters of
   the last resolution level. */
  itkSetMacro( InitialTransformParametersOfNextLevel, ParametersType );
  itkGetConstReferenceMacro( InitialTransformParametersOfNextLevel, ParametersType );

  /** Method to return the latest modified time of this object or
   * any of its cached ivars */
  virtual unsigned long GetMTime() const;

  /** Initialize by setting the interconnects between the components.
      This method is executed at every level of the pyramid with the
      values corresponding to this resolution
   */
  void Initialize() throw (ExceptionObject);

protected:
  MultiResolutionMultiImageToImageRegistrationMethod();
  virtual ~MultiResolutionMultiImageToImageRegistrationMethod() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  void  GenerateData ();

  /** Compute the size of the fixed region for each level of the pyramid. */
  void PreparePyramids( void );

  /** Set the current level to be processed */
  itkSetMacro( CurrentLevel, unsigned long );

private:
  MultiResolutionMultiImageToImageRegistrationMethod(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  FixedMultiImageType              m_CurrentFixedMultiImage;

  MovingImagePyramidPointer        m_MovingImagePyramid;
  FixedMultiImagePyramidType       m_FixedMultiImagePyramid;

  ParametersType                   m_InitialTransformParametersOfNextLevel;

  FixedMultiImageRegionPyramidType m_FixedMultiImageRegionPyramid;
  FixedMultiImageRegionType        m_CurrentFixedMultiImageRegion;

  unsigned long                    m_NumberOfLevels;
  unsigned long                    m_CurrentLevel;

  bool                             m_Stop;

  ScheduleType                     m_FixedImagePyramidSchedule;
  ScheduleType                     m_MovingImagePyramidSchedule;

  bool                             m_ScheduleSpecified;
  bool                             m_NumberOfLevelsSpecified;

};


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiResolutionMultiImageToImageRegistrationMethod.txx"
#endif

#endif
