#ifndef __itkGradientDifferenceMultiImageToImageMetric_h
#define __itkGradientDifferenceMultiImageToImageMetric_h

#include "itkMultiImageToImageMetric.h"
#include "itkGradientDifferenceSingleImageToImageMetric.h"

namespace itk{

/** \class GradientDifferenceMultiImageToImageMetric
 * \brief Version of the Gradient Difference metric for registrations with multiple fixed images.
 *
 */
  template < class TFixedImage, class TMovingImage >
  class ITK_EXPORT GradientDifferenceMultiImageToImageMetric :
  public MultiImageToImageMetric<TFixedImage,TMovingImage>
  {
  public:
/** Standard class typedefs */
    typedef GradientDifferenceMultiImageToImageMetric         Self;
    typedef MultiImageToImageMetric<TFixedImage,TMovingImage> Superclass;
    typedef SmartPointer<Self>                                Pointer;
    typedef SmartPointer<const Self>                          ConstPointer;

/** Method for creation through the object factory. */
    itkNewMacro(Self);

/** Run-time type information (and related methods). */
    itkTypeMacro(GradientDifferenceMultiImageToImageMetric, MultiImageToImageMetric);

/** Types defined by this class */
    typedef GradientDifferenceSingleImageToImageMetric< TFixedImage, TMovingImage > MetricType;

  protected:
    GradientDifferenceMultiImageToImageMetric() { itkDebugMacro("Constructor"); }
    virtual ~GradientDifferenceMultiImageToImageMetric() {};

    itkNewSingleMetricMacro(MetricType);

  private:
    GradientDifferenceMultiImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&);          //purposely not implemented
  };
} //namespace itk

#endif