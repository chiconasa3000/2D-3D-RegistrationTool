#ifndef __itkMeanSquaresMultiImageToImageMetric_h
#define __itkMeanSquaresMultiImageToImageMetric_h

#include "itkMultiImageToImageMetric.h"
#include <itkMeanSquaresImageToImageMetric.h>

namespace itk{

/** \class MeanSquaresMultiImageToImageMetric
 * \brief Version of the mean squares metric for registrations with multiple fixed images.
 * 
*/
  template < class TFixedImage, class TMovingImage >
  class ITK_EXPORT MeanSquaresMultiImageToImageMetric :
  public MultiImageToImageMetric<TFixedImage,TMovingImage>
  {
  public:
/** Standard class typedefs */
    typedef MeanSquaresMultiImageToImageMetric                Self;
    typedef MultiImageToImageMetric<TFixedImage,TMovingImage> Superclass;
    typedef SmartPointer<Self>                                Pointer;
    typedef SmartPointer<const Self>                          ConstPointer;

/** Method for creation through the object factory. */
    itkNewMacro(Self);

/** Run-time type information (and related methods). */
    itkTypeMacro(MeanSquaresMultiImageToImageMetric, MultiImageToImageMetric);

/** Types defined by this class */
    typedef MeanSquaresImageToImageMetric<TFixedImage,TMovingImage> MetricType;

  protected:
    MeanSquaresMultiImageToImageMetric() { itkDebugMacro("Constructor"); }
    virtual ~MeanSquaresMultiImageToImageMetric() {};

    itkNewSingleMetricMacro(MetricType);

  private:
    MeanSquaresMultiImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&);          //purposely not implemented
  };
} //namespace itk

#endif