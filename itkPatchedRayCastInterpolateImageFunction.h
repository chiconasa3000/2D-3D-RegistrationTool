#ifndef __itkPatchedRayCastInterpolateImageFunction_h
#define __itkPatchedRayCastInterpolateImageFunction_h

#include <itkInterpolateImageFunction.h>
//#include <itkEuler3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkTransform.h>
#include <itkVector.h>

namespace itk
{

/** \class PatchedRayCastInterpolateImageFunction
 * \brief Projective interpolation of an image at specified positions.
 *
 * PatchedRayCastInterpolateImageFunction casts rays through a 3-dimensional
 * image and uses bilinear interpolation to integrate each plane of
 * voxels traversed.
 *
 * This class is almost identical to the ITK's original
 * RayCastInterpolateImageFunction, but was patched to support more complex
 * geometric definitions of the fixed and moving images. Some of the applied
 * fixes were proposed by Jian Wu and Aviv Hurvitz on the following posts:
 *
 * http://public.kitware.com/pipermail/insight-users/2006-March/017265.html
 *
 * http://blog.gmane.org/gmane.comp.lib.itk.user/month=20080701
 *
 * In addition, an extra patch was applied on the Evaluate method to handle
 * the input image's orientation correctly, something which was not done by
 * the original class.
 *
 * \warning This interpolator works for 3-dimensional images only.
 *
 */
template <class TInputImage, class TCoordRep = double>
class ITK_EXPORT PatchedRayCastInterpolateImageFunction :
    public InterpolateImageFunction<TInputImage,TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef PatchedRayCastInterpolateImageFunction                 Self;
  typedef InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Constants for the image dimensions */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /**
   * Type of the Transform Base class
   * The fixed image should be a 3D image
   */
  typedef Transform<TCoordRep,3,3> TransformType;

  typedef typename TransformType::Pointer            TransformPointer;
  typedef typename TransformType::InputPointType     InputPointType;
  typedef typename TransformType::OutputPointType    OutputPointType;
  typedef typename TransformType::ParametersType     TransformParametersType;
  typedef typename TransformType::JacobianType       TransformJacobianType;

  typedef typename Superclass::InputPixelType        PixelType;

  typedef typename TInputImage::SizeType             SizeType;

  typedef Vector<TCoordRep, 3>                       DirectionType;

  /**  Type of the Interpolator Base class */
  typedef InterpolateImageFunction<TInputImage,TCoordRep> InterpolatorType;

  typedef typename InterpolatorType::Pointer         InterpolatorPointer;


  /** Run-time type information (and related methods). */
  itkTypeMacro(PatchedRayCastInterpolateImageFunction, InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Post-transform typedef */
  typedef typename itk::Similarity3DTransform<TCoordRep> PostTransformType;
  typedef typename PostTransformType::Pointer       PostTransformPointer;

  /** \brief
   * Interpolate the image at a point position.
   *
   * Returns the interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method.
   */
  virtual OutputType Evaluate( const PointType& point ) const;

  /** Interpolate the image at a continuous index position
   *
   * Returns the interpolated image intensity at a
   * specified index position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * Subclasses must override this method.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method.
   */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType &index ) const;


  /** Connect the Transform. */
  itkSetObjectMacro( Transform, TransformType );
  /** Get a pointer to the Transform.  */
  itkGetObjectMacro( Transform, TransformType );

  /** Connect the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );
  /** Get a pointer to the Interpolator.  */
  itkGetObjectMacro( Interpolator, InterpolatorType );

  /** Connect the Interpolator. */
  itkSetMacro( FocalPoint, InputPointType );
  /** Get a pointer to the Interpolator.  */
  itkGetConstMacro( FocalPoint, InputPointType );

  /** Connect the Transform. */
  itkSetMacro( Threshold, double );
  /** Get a pointer to the Transform.  */
  itkGetConstMacro( Threshold, double );

  /** Connect the Image */
  virtual void SetInputImage( const InputImageType* );

  /** Check if a point is inside the image buffer.
   * \warning For efficiency, no validity checking of
   * the input image pointer is done. */
  inline bool IsInsideBuffer( const PointType & ) const
    {
    return true;
    }
  bool IsInsideBuffer( const ContinuousIndexType &  ) const
    {
    return true;
    }
  bool IsInsideBuffer( const IndexType &  ) const
    {
    return true;
    }

protected:

  /// Constructor
  PatchedRayCastInterpolateImageFunction();

  /// Destructor
  ~PatchedRayCastInterpolateImageFunction(){};

  /// Print the object
  void PrintSelf(std::ostream& os, Indent indent) const;

  /// Transformation used to calculate the new focal point position
  TransformPointer m_Transform;

  /// The focal point or position of the ray source
  InputPointType m_FocalPoint;

  /// The threshold above which voxels along the ray path are integrated.
  double m_Threshold;

  /// Pointer to the interpolator
  InterpolatorPointer m_Interpolator;


private:
  PatchedRayCastInterpolateImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  PostTransformPointer m_PostTransform;

};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPatchedRayCastInterpolateImageFunction.txx"
#endif

#endif
