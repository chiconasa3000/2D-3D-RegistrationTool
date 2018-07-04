#ifndef __itkNormalizedGradientCorrelationMultiImageToImageMetric_h
#define __itkNormalizedGradientCorrelationMultiImageToImageMetric_h

#include "itkMultiImageToImageMetric.h"
#include "itkNormalizedGradientCorrelationImageToImageMetric.h"

namespace itk{

/** \class NormalizedGradientCorrelationMultiImageToImageMetric
 *
 * \brief Version of the normalized gradient correlation metric for registrations with multiple fixed images.
 */
  template < class TFixedImage, class TMovingImage >
  class ITK_EXPORT NormalizedGradientCorrelationMultiImageToImageMetric :
  public MultiImageToImageMetric<TFixedImage,TMovingImage>
  {
  public:
/** Standard class typedefs */
    typedef NormalizedGradientCorrelationMultiImageToImageMetric                Self;
    typedef MultiImageToImageMetric<TFixedImage,TMovingImage> Superclass;
    typedef SmartPointer<Self>                                Pointer;
    typedef SmartPointer<const Self>                          ConstPointer;

/** Method for creation through the object factory. */
    itkNewMacro(Self);

/** Run-time type information (and related methods). */
    itkTypeMacro(NormalizedGradientCorrelationMultiImageToImageMetric, MultiImageToImageMetric);

/** Types defined by this class */
    typedef NormalizedGradientCorrelationImageToImageMetric< TFixedImage, TMovingImage > MetricType;

  protected:
    NormalizedGradientCorrelationMultiImageToImageMetric() {  itkDebugMacro("Constructor"); }
    virtual ~NormalizedGradientCorrelationMultiImageToImageMetric() {};

    itkNewSingleMetricMacro(MetricType);

  private:
    NormalizedGradientCorrelationMultiImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&);          //purposely not implemented
  };
} //namespace itk

#endif