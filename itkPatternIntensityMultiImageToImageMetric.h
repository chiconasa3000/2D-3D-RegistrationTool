#ifndef __itkPatternIntensityMultiImageToImageMetric_h
#define __itkPatternIntensityMultiImageToImageMetric_h

#include "itkMultiImageToImageMetric.h"
#include "itkPatternIntensityImageToImageMetric.h"

namespace itk{

/** \class PatternIntensityMultiImageToImageMetric
 *
 * \brief Version of the pattern intensity metric for registrations with multiple fixed images.
 *
 * This class implements the Pattern Intensity metric for registrations with
 * multiple fixed images. The class includes the SetRadius and SetSigma
 * methods, which sets the Radius and Sigma parameters on each of the
 * individual metrics contained within the class (please see the
 * PatternIntensityImageToImageMetric class for a brief description of these
 * parameters).
 *
 */
  template < class TFixedImage, class TMovingImage >
  class ITK_EXPORT PatternIntensityMultiImageToImageMetric :
  public MultiImageToImageMetric<TFixedImage,TMovingImage>
  {
  public:
/** Standard class typedefs */
    typedef PatternIntensityMultiImageToImageMetric           Self;
    typedef MultiImageToImageMetric<TFixedImage,TMovingImage> Superclass;
    typedef SmartPointer<Self>                                Pointer;
    typedef SmartPointer<const Self>                          ConstPointer;

/** Method for creation through the object factory. */
    itkNewMacro(Self);

/** Run-time type information (and related methods). */
    itkTypeMacro(PatternIntensityMultiImageToImageMetric, MultiImageToImageMetric);

/** Types defined by this class */
    typedef PatternIntensityImageToImageMetric< TFixedImage, TMovingImage > MetricType;
    typedef typename MetricType::Pointer                                    MetricPointer;

/** Set or Get the radius */
    itkSetMacro(Radius,int);
    itkGetConstReferenceMacro(Radius,int);

/** Set or Get the sigma constant */
    itkSetMacro(Sigma,double);
    itkGetConstReferenceMacro(Sigma,double);

  protected:
    PatternIntensityMultiImageToImageMetric();
    virtual ~PatternIntensityMultiImageToImageMetric() {};
    void PrintSelf(std::ostream& os, Indent indent) const;

    //Cannot use the itkNewSingleMetric macro, as new single metrics need
    //proper initializations of NoiseConstant and NeighborhoodRadius
    typedef typename Superclass::MetricType  BaseMetricType;
    typedef typename BaseMetricType::Pointer BaseMetricPointer;
    virtual BaseMetricPointer NewSingleMetric();

    virtual void DoNumberRevision(void) const throw ( ExceptionObject );

  private:
    PatternIntensityMultiImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&);                          //purposely not implemented

    int m_Radius;
    double m_Sigma;
  };
} //namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPatternIntensityMultiImageToImageMetric.txx"
#endif

#endif
