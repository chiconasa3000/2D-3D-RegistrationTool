/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineDeformableTransform.h,v $
  Language:  C++
  Date:      $Date: 2006/11/03 20:09:08 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __BSplineDeformableTransformOpt_h
#define __BSplineDeformableTransformOpt_h

#include <iostream>
//#include "itkTransform.h"
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkBSplineInterpolationWeightFunction.h"

namespace itk
{

/** \class BSplineDeformableTransform
   * \brief Bspline with optimization
   *
   * This class computes the jacobian locally.
   * It does not need to store the jacobian field.
   * 
   * GetJacobian( point, indexes, weights )
   * returns the nonzero indexes and coressponding weights 
   * for the jacobian at a given point.
   *
   * This class is backward compatible as it also provides
   * GetJacobian( point)
   *
   * \ingroup Transforms
 */

template <
    class TScalarType = double,          // Data type for scalars
    unsigned int NDimensions = 3,        // Number of dimensions
    unsigned int VSplineOrder = 3 >      // Spline order
class ITK_EXPORT BSplineDeformableTransformOpt :
          public Transform< TScalarType, NDimensions, NDimensions >
{
public:
  /** Standard class typedefs. */
  typedef BSplineDeformableTransformOpt                         Self;
  typedef Transform< TScalarType, NDimensions, NDimensions > Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;
      
  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( BSplineDeformableTransformOpt, Transform );

  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);

  /** The BSpline order. */
  itkStaticConstMacro(SplineOrder, unsigned int, VSplineOrder);

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard parameters container. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Standard Jacobian container. */
  typedef typename Superclass::JacobianType JacobianType;

  /** Standard vector type for this class. */
  typedef Vector<TScalarType,  itkGetStaticConstMacro(SpaceDimension)> InputVectorType;
  typedef Vector<TScalarType,  itkGetStaticConstMacro(SpaceDimension)> OutputVectorType;

  /** Standard covariant vector type for this class. */
  typedef CovariantVector<TScalarType,  itkGetStaticConstMacro(SpaceDimension)> InputCovariantVectorType;
  typedef CovariantVector<TScalarType,  itkGetStaticConstMacro(SpaceDimension)> OutputCovariantVectorType;
  
  /** Standard vnl_vector type for this class. */
  typedef vnl_vector_fixed<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputVnlVectorType;
  typedef vnl_vector_fixed<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputVnlVectorType;
  
  /** Standard coordinate point type for this class. */
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputPointType;
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputPointType;
  
  /** This method sets the parameters of the transform.
   * For a BSpline deformation transform, the parameters are the BSpline 
   * coefficients on a sparse grid. 
   * 
   * The parameters are N number of N-D grid of coefficients. Each N-D grid 
   * is represented as a flat array of doubles 
   * (in the same configuration as an itk::Image).
   * The N arrays are then concatenated to form one parameter array.
   *
   * For efficiency, this transform does not make a copy of the parameters.
   * It only keeps a pointer to the input parameters. It assumes that the memory
   * is managed by the caller. Use SetParametersByValue to force the transform
   * to call copy the parameters.
   *
   * This method wraps each grid as itk::Image's using the user specified
   * grid region, spacing and origin.
   * NOTE: The grid region, spacing and origin must be set first.
   *
   */

  //----------------  implemented by bsplinebasetransform -----------------//
  void SetParameters(const ParametersType & parameters);
  
  /** This method sets the fixed parameters of the transform.
   * For a BSpline deformation transform, the parameters are the following:
   *    Grid Size, Grid Origin, and Grid Spacing
   * 
   * The fixed parameters are the three times the size of the templated 
   * dimensions.
   * This function has the effect of make the following calls:
   *       transform->SetGridSpacing( spacing );
   *       transform->SetGridOrigin( origin );
   *       transform->SetGridRegion( bsplineRegion );
   *
   * This function was added to allow the transform to work with the 
   * itkTransformReader/Writer I/O filters.
   *
   */

  //----------------  implemented by bsplinebasetransform -----------------//
  void SetFixedParameters(const ParametersType & parameters);

  /** This method sets the parameters of the transform.
   * For a BSpline deformation transform, the parameters are the BSpline 
   * coefficients on a sparse grid. 
   * 
   * The parameters are N number of N-D grid of coefficients. Each N-D grid 
   * is represented as a flat array of doubles 
   * (in the same configuration as an itk::Image).
   * The N arrays are then concatenated to form one parameter array.
   *
   * This methods makes a copy of the parameters while for
   * efficiency the SetParameters method does not.
   *
   * This method wraps each grid as itk::Image's using the user specified
   * grid region, spacing and origin.
   * NOTE: The grid region, spacing and origin must be set first.
   *
   */

  //----------------  implemented by bsplinebasetransform -----------------//
  void SetParametersByValue(const ParametersType & parameters);

  /** This method can ONLY be invoked AFTER calling SetParameters(). 
   *  This restriction is due to the fact that the BSplineDeformableTransform
   *  does not copy the array of paramters internally, instead it keeps a 
   *  pointer to the user-provided array of parameters. This method is also
   *  in violation of the const-correctness of the parameters since the 
   *  parameter array has been passed to the transform on a 'const' basis but
   *  the values get modified when the user invokes SetIdentity().
   */
  //----------------  implemented by bsplinebasetransform -----------------//
  void SetIdentity();

  /** Get the Transformation Parameters. */
  //----------------  implemented by bsplinebasetransform -----------------//
  virtual const ParametersType& GetParameters(void) const;

  /** Get the Transformation Fixed Parameters. */
  //----------------  implemented by bsplinebasetransform -----------------//
  virtual const ParametersType& GetFixedParameters(void) const;
  
  /** Parameters as SpaceDimension number of images. */
  typedef typename ParametersType::ValueType                      PixelType;
  typedef Image<PixelType,itkGetStaticConstMacro(SpaceDimension)> ImageType;
  typedef typename ImageType::Pointer                             ImagePointer;

  /** Get the array of coefficient images. */

  //----------------  implemented by bsplinebasetransform -----------------//
  virtual ImagePointer * GetCoefficientImage()
    { return m_CoefficientImage; }

  /** Set the array of coefficient images.
   *
   * This is an alternative API for setting the BSpline coefficients
   * as an array of SpaceDimension images. The grid region spacing 
   * and origin is taken from the first image. It is assume that
   * the buffered region of all the subsequent images are the same 
   * as the first image. Note that no error checking is done.
   *
   * Warning: use either the SetParameters() or SetCoefficientImage()
   * API. Mixing the two modes may results in unexpected results.
   *
   */
   //----------------  implemented by bsplinebasetransform -----------------//
  virtual void SetCoefficientImage( ImagePointer images[] );  

  /** Typedefs for specifying the extend to the grid. */
  typedef ImageRegion<itkGetStaticConstMacro(SpaceDimension)>    RegionType;
  
  typedef typename RegionType::IndexType  IndexType;
  typedef typename RegionType::SizeType   SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::PointType   OriginType;

  /** This method specifies the region over which the grid resides. */
  virtual void SetGridRegion( const RegionType& region );
  itkGetMacro( GridRegion, RegionType );
  itkGetConstMacro( GridRegion, RegionType );

  /** This method specifies the grid spacing or resolution. */
  virtual void SetGridSpacing( const SpacingType& spacing );
  itkGetMacro( GridSpacing, SpacingType );
  itkGetConstMacro( GridSpacing, SpacingType );

  /** This method specifies the grid origin. */
  virtual void SetGridOrigin( const OriginType& origin );
  itkGetMacro( GridOrigin, OriginType );
  itkGetConstMacro( GridOrigin, OriginType );

  /** Typedef of the bulk transform. */
  typedef Transform<ScalarType, itkGetStaticConstMacro(SpaceDimension), itkGetStaticConstMacro(SpaceDimension)> BulkTransformType;
  typedef typename BulkTransformType::ConstPointer   BulkTransformPointer;

  /** This method specifies the bulk transform to be applied. 
   * The default is the identity transform.
   */
  itkSetConstObjectMacro( BulkTransform, BulkTransformType );
  itkGetConstObjectMacro( BulkTransform, BulkTransformType );

  /** Transform points by a BSpline deformable transformation. */
  //esta funcion proviene de bsplinebasetransform asi que se obvia
  //pero puede considerarse la implementacion de aqui OJO
  OutputPointType  TransformPoint(const InputPointType  &point ) const;

  /** Interpolation weights function type. */
  typedef BSplineInterpolationWeightFunction<ScalarType,itkGetStaticConstMacro(SpaceDimension), itkGetStaticConstMacro(SplineOrder)> WeightsFunctionType;
  typedef typename WeightsFunctionType::WeightsType       WeightsType;
  typedef typename WeightsFunctionType::ContinuousIndexType ContinuousIndexType;

  /** Parameter index array type. */
  typedef Array<unsigned long> ParameterIndexArrayType;

  /** Transform points by a BSpline deformable transformation. 
   * On return, weights contains the interpolation weights used to compute the 
   * deformation and indices of the x (zeroth) dimension coefficient parameters
   * in the support region used to compute the deformation.
   * Parameter indices for the i-th dimension can be obtained by adding
   * ( i * this->GetNumberOfParametersPerDimension() ) to the indices array.
   */
  //esta funcion ya ha sido implementada en la actual clase que es del padre OJO
  //pero con dos parametros adicionales   ParameterIndexArrayType &   indices y un bool &    inside 
  virtual void TransformPoint( const InputPointType & inputPoint, OutputPointType & outputPoint, WeightsType & weights ) const;

  /** Get number of weights. */
  unsigned long GetNumberOfWeights() const
    { return m_WeightsFunction->GetNumberOfWeights(); }

  /** Method to transform a vector - 
   *  not applicable for this type of transform. */
  virtual OutputVectorType TransformVector(const InputVectorType &) const
    { 
    itkExceptionMacro(<< "Method not applicable for deformable transform." );
    return OutputVectorType(); 
    }

  /** Method to transform a vnl_vector - 
   *  not applicable for this type of transform */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const
    { 
    itkExceptionMacro(<< "Method not applicable for deformable transform. ");
    return OutputVnlVectorType(); 
    }

  /** Method to transform a CovariantVector - 
   *  not applicable for this type of transform */
 virtual OutputCovariantVectorType TransformCovariantVector(
    const InputCovariantVectorType &) const
    { 
    itkExceptionMacro(<< "Method not applicable for deformable transfrom. ");
    return OutputCovariantVectorType(); 
    } 
   //This Functions must be used but the current implementation
 //used his own getjacobian so the superclass transform
 //needs computejacobian functions
 //so I declared this as virtual function
 //because we don't use this functions.
 virtual void ComputeJacobianWithRespectToParameters(const InputPointType &, JacobianType &) const
 {
    itkExceptionMacro(<< "Method for version 4 in ITK ");
    return; 
 }
virtual void ComputeJacobianWithRespectToPosition(const InputPointType & x, JacobianType &jacobian) const
 {
    itkExceptionMacro(<< "Method for version 4 in ITK ");
    return; 
 }
 
  /** Compute the Jacobian Matrix of the transformation at one point */
  //--------- accedido mediante void ComputeJacobianWithRespectToParameters (const InputPointType &, JacobianType &) 
  const JacobianType& GetJacobian(const InputPointType  &point ) const;
  
  //-------- accedido mediante void ComputeJacobianFromBSplineWeightsWithRespectToPosition (const InputPointType &, WeightsType &, ParameterIndexArrayType &) const 
  void GetJacobian( const InputPointType & point, ParameterIndexArrayType& indexes, WeightsType & weights );
  
  /** Return the number of parameters that completely define the Transfom */
   //----------------  implemented by bsplinebasetransform -----------------//
  virtual unsigned long GetNumberOfParameters(void) const;

  /** Return the number of parameters per dimension */
   //----------------  implemented by bsplinebasetransform -----------------// 
  //OJO podemos usar su propia implementacion
  unsigned long GetNumberOfParametersPerDimension(void) const;

  /** Return the region of the grid wholly within the support region */
  itkGetConstReferenceMacro( ValidRegion, RegionType );

  /** Indicates that this transform is linear. That is, given two
   * points P and Q, and scalar coefficients a and b, then
   *
   *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
   */
  virtual bool IsLinear() const { return false; }

  /** Get number of weights which affect a points value at a particular place */
   //----------------  implemented by bsplinebasetransform -----------------// 
  //OJO podemos usar su propia implementacion
  unsigned int GetNumberOfAffectedWeights() const;

protected:
  /** Print contents of an BSplineDeformableTransform. */
  void PrintSelf(std::ostream &os, Indent indent) const;


  BSplineDeformableTransformOpt();
  virtual ~BSplineDeformableTransformOpt();

  /** Allow subclasses to access and manipulate the weights function. */
  itkSetObjectMacro( WeightsFunction, WeightsFunctionType );
  itkGetObjectMacro( WeightsFunction, WeightsFunctionType );

  /** Wrap flat array into images of coefficients. */
  //----------------  implemented by bsplinebasetransform -----------------// 
  //OJO podemos usar su propia implementacion
  void WrapAsImages();

private:
  BSplineDeformableTransformOpt(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The bulk transform. */
  BulkTransformPointer  m_BulkTransform;

  /** Variables defining the coefficient grid extend. */
  RegionType    m_GridRegion;
  SpacingType   m_GridSpacing;
  OriginType    m_GridOrigin;
  
  RegionType    m_ValidRegion;

  /** Variables defining the interpolation support region. */
  unsigned long m_Offset;
  bool          m_SplineOrderOdd;
  SizeType      m_SupportSize;
  IndexType     m_ValidRegionLast;
  
  /** Array holding images wrapped from the flat parameters. */
  ImagePointer   m_WrappedImage[NDimensions];

/*His own jacobian variable*/
JacobianType m_Jacobian;

  /** Array of images representing the B-spline coefficients 
   *  in each dimension. */
  ImagePointer   m_CoefficientImage[NDimensions];

  /** Jacobian as SpaceDimension number of images. */
  typedef typename JacobianType::ValueType JacobianPixelType;
  typedef Image<JacobianPixelType, itkGetStaticConstMacro(SpaceDimension)> JacobianImageType;

  typename JacobianImageType::Pointer m_JacobianImage[NDimensions];

  /** Keep track of last support region used in computing the Jacobian
   * for fast resetting of Jacobian to zero.
   */
  mutable IndexType m_LastJacobianIndex;

  /** Keep a pointer to the input parameters. */
  const ParametersType *  m_InputParametersPointer;

  /** Internal parameters buffer. */
  ParametersType          m_InternalParametersBuffer;

  /** Pointer to function used to compute Bspline interpolation weights. */
  typename WeightsFunctionType::Pointer  m_WeightsFunction;

  /** Check if a continuous index is inside the valid region. */
  //----------------  implemented by bsplinebasetransform -----------------// 
  //OJO podemos usar su propia implementacion
  bool InsideValidRegion( const ContinuousIndexType& index ) const;


}; //class BSplineDeformableTransformOpt


}  // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "BSplineDeformableTransformOpt.txx"
#endif

#endif /* __itkBSplineDeformableTransform_h */
