/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineDeformableTransform.txx,v $
  Language:  C++
  Date:      $Date: 2006/07/11 15:13:26 $
  Version:   $Revision: 1.29 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __BSplineDeformableTransformOpt_txx
#define __BSplineDeformableTransformOpt_txx

#include "BSplineDeformableTransformOpt.h"
#include "itkContinuousIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkIdentityTransform.h"

namespace itk
{

// Constructor with default arguments
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::BSplineDeformableTransformOpt():Superclass(0)
{

  // Instantiate a weights function
  m_WeightsFunction = WeightsFunctionType::New();
  m_SupportSize = m_WeightsFunction->GetSupportSize();

  // Instantiate an identity transform
  typedef IdentityTransform<ScalarType,SpaceDimension> IdentityTransformType;
  typename IdentityTransformType::Pointer id = IdentityTransformType::New();
  m_BulkTransform = id;

  // Default grid size is zero
  typename RegionType::SizeType size;
  typename RegionType::IndexType index;
  size.Fill( 0 );
  index.Fill( 0 );
  m_GridRegion.SetSize( size );
  m_GridRegion.SetIndex( index );

  m_GridOrigin.Fill( 0.0 );  // default origin is all zeros
  m_GridSpacing.Fill( 1.0 ); // default spacing is all ones

  m_InputParametersPointer = NULL;
  m_InternalParametersBuffer = ParametersType(0);

  // Initialize coeffient images
  for ( unsigned int j = 0; j < SpaceDimension; j++ )
    {
    m_WrappedImage[j] = ImageType::New();
    m_WrappedImage[j]->SetRegions( m_GridRegion );
    m_WrappedImage[j]->SetOrigin( m_GridOrigin.GetDataPointer() );
    m_WrappedImage[j]->SetSpacing( m_GridSpacing.GetDataPointer() );
    m_CoefficientImage[j] = NULL;
    }

  // Setup variables for computing interpolation
  m_Offset = SplineOrder / 2;
  if ( SplineOrder % 2 ) 
    {
    m_SplineOrderOdd = true;
    }
  else
    {
    m_SplineOrderOdd = false;
    }
  m_ValidRegion = m_GridRegion;

  // Initialize jacobian images
  for ( unsigned int j = 0; j < SpaceDimension; j++ )
    {
    m_JacobianImage[j] = ImageType::New();
    m_JacobianImage[j]->SetRegions( m_GridRegion );
    m_JacobianImage[j]->SetOrigin( m_GridOrigin.GetDataPointer() );
    m_JacobianImage[j]->SetSpacing( m_GridSpacing.GetDataPointer() );
    }

  /** Fixed Parameters store the following information:
   *     Grid Size
   *     Grid Origin
   *     Grid Spacing
   *  The size of these is equal to the  NInputDimensions */
  this->m_FixedParameters.SetSize ( NDimensions * 3 );
  this->m_FixedParameters.Fill ( 0.0 );
  
  m_LastJacobianIndex = m_ValidRegion.GetIndex();
  
}
    

// Destructor
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::~BSplineDeformableTransformOpt()
{

}


// Get the number of parameters
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
unsigned long
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::GetNumberOfParameters(void) const
{

  // The number of parameters equal SpaceDimension * number of
  // of pixels in the grid region.
  return ( static_cast<unsigned long>( SpaceDimension ) * static_cast<unsigned long>( m_GridRegion.GetNumberOfPixels() ) );

}


// Get the number of parameters per dimension
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
unsigned long
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::GetNumberOfParametersPerDimension(void) const
{
  // The number of parameters per dimension equal number of
  // of pixels in the grid region.
  return ( static_cast<unsigned long>( m_GridRegion.GetNumberOfPixels() ) );

}


// Set the grid region
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetGridRegion( const RegionType& region )
{
  if ( m_GridRegion != region )
    {

    m_GridRegion = region;

    // set regions for each coefficient and jacobian image
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_WrappedImage[j]->SetRegions( m_GridRegion );
      m_JacobianImage[j]->SetRegions( m_GridRegion );
      }

    // Set the valid region
    // If the grid spans the interval [start,last].
    // The valid interval for evaluation is [start+offset,last-offset]
    // when spline order is even.
    // The valid interval for evaluation is [start+offset,last-offset)
    // when spline order is odd.
    // Where offset = vcl_floor(spline / 2 ).
    // Note that the last pixel is not included in the valid region
    // with odd spline orders.
    typename RegionType::SizeType size = m_GridRegion.GetSize();
    typename RegionType::IndexType index = m_GridRegion.GetIndex();
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      index[j] += 
        static_cast< typename RegionType::IndexValueType >( m_Offset );
      size[j] -= 
        static_cast< typename RegionType::SizeValueType> ( 2 * m_Offset );
      m_ValidRegionLast[j] = index[j] +
        static_cast< typename RegionType::IndexValueType >( size[j] ) - 1;
      }
    m_ValidRegion.SetSize( size );
    m_ValidRegion.SetIndex( index );

    this->Modified();
    }
}



// Set the grid spacing
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetGridSpacing( const SpacingType& spacing )
{
  if ( m_GridSpacing != spacing )
    {
    m_GridSpacing = spacing;

    // set spacing for each coefficient and jacobian image
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_WrappedImage[j]->SetSpacing( m_GridSpacing.GetDataPointer() );
      m_JacobianImage[j]->SetSpacing( m_GridSpacing.GetDataPointer() );
      }

    this->Modified();
    }

}



// Set the grid origin
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetGridOrigin( const OriginType& origin )
{
  if ( m_GridOrigin != origin )
    {
    m_GridOrigin = origin;

    // set spacing for each coefficient and jacobianimage
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_WrappedImage[j]->SetOrigin( m_GridOrigin.GetDataPointer() );
      m_JacobianImage[j]->SetOrigin( m_GridOrigin.GetDataPointer() );
      }

    this->Modified();
    }

}


// Set the parameters
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetIdentity()
{
  if( m_InputParametersPointer )
    {
    ParametersType * parameters =
      const_cast<ParametersType *>( m_InputParametersPointer );
    parameters->Fill( 0.0 );
    this->Modified();
    }
  else 
    {
    itkExceptionMacro( << "Input parameters for the spline haven't been set ! "
       << "Set them using the SetParameters or SetCoefficientImage method first." );
    }
}

// Set the parameters

template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetParameters( const ParametersType & parameters )
{

  // check if the number of parameters match the
  // expected number of parameters
  if ( parameters.Size() != this->GetNumberOfParameters() )
    {
    itkExceptionMacro(<<"Mismatched between parameters size "
                      << parameters.size() 
                      << " and region size "
                      << m_GridRegion.GetNumberOfPixels() );
    }

  // Clean up buffered parameters
  m_InternalParametersBuffer = ParametersType( 0 );

  // Keep a reference to the input parameters
  m_InputParametersPointer = &parameters;

  // Wrap flat array as images of coefficients
  this->WrapAsImages();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();
}

// Set the Fixed Parameters
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetFixedParameters( const ParametersType & parameters )
{
 
  // check if the number of parameters match the
  // expected number of parameters
  if ( parameters.Size() != NDimensions*3 )
    {
    itkExceptionMacro(<< "Mismatched between parameters size "
                      << parameters.size() 
                      << " and number of fixed parameters "
                      << NDimensions*3 );
    }

  /********************************************************** 
    Fixed Parameters store the following information:
        Grid Size
        Grid Origin
        Grid Spacing
     The size of these is equal to the  NInputDimensions
  **********************************************************/
  
  /*** Set the Grid Parameters ***/
  SizeType   gridSize;
  for (unsigned int i=0;i<NDimensions;i++)
    {
    gridSize[i] = static_cast<int> (parameters[i]);
    }
  RegionType bsplineRegion;
  bsplineRegion.SetSize( gridSize );
  
  /*** Set the Origin Parameters ***/
  OriginType origin;
  for (unsigned int i=0;i<NDimensions;i++)
    {
    origin[i] = parameters[NDimensions+i];
    }
  
  /*** Set the Spacing Parameters ***/
  SpacingType spacing;
  for (unsigned int i=0;i<NDimensions;i++)
    {
    spacing[i] = parameters[2*NDimensions+i];
    }

  
  this->SetGridSpacing( spacing );
  this->SetGridOrigin( origin );
  this->SetGridRegion( bsplineRegion );

  this->Modified();
}


// Wrap flat parameters as images
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::WrapAsImages()
{

  /**
   * Wrap flat parameters array into SpaceDimension number of ITK images
   * NOTE: For efficiency, parameters are not copied locally. The parameters
   * are assumed to be maintained by the caller.
   */
  PixelType * dataPointer =
    const_cast<PixelType *>(( m_InputParametersPointer->data_block() ));
  unsigned int numberOfPixels = m_GridRegion.GetNumberOfPixels();

  for ( unsigned int j = 0; j < SpaceDimension; j++ )
    {
    m_WrappedImage[j]->GetPixelContainer()->
      SetImportPointer( dataPointer, numberOfPixels );
    dataPointer += numberOfPixels;
    m_CoefficientImage[j] = m_WrappedImage[j];
    }

    
  /**
     * Allocate memory for Jacobian and wrap into SpaceDimension number
     * of ITK images
   */
    this->m_Jacobian.set_size( SpaceDimension, this->GetNumberOfParameters() );
    this->m_Jacobian.Fill( NumericTraits<JacobianPixelType>::Zero );
    m_LastJacobianIndex = m_ValidRegion.GetIndex();
    JacobianPixelType * jacobianDataPointer = this->m_Jacobian.data_block();

    for ( unsigned int j = 0; j < SpaceDimension; j++ )
    {
      m_JacobianImage[j]->GetPixelContainer()->
          SetImportPointer( jacobianDataPointer, numberOfPixels );
      jacobianDataPointer += this->GetNumberOfParameters() + numberOfPixels;
    }
    
}


// Set the parameters by value
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetParametersByValue( const ParametersType & parameters )
{

  // check if the number of parameters match the
  // expected number of parameters
  if ( parameters.Size() != this->GetNumberOfParameters() )
    {
    itkExceptionMacro(<<"Mismatched between parameters size "
                      << parameters.size() 
                      << " and region size "
                      << m_GridRegion.GetNumberOfPixels() );
    }

  // copy it
  m_InternalParametersBuffer = parameters;
  m_InputParametersPointer = &m_InternalParametersBuffer;

  // wrap flat array as images of coefficients
  this->WrapAsImages();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

}



// Get the parameters
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
const 
typename BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::ParametersType &
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::GetParameters( void ) const
{
  /** NOTE: For efficiency, this class does not keep a copy of the parameters - 
   * it just keeps pointer to input parameters. 
   */
  return (*m_InputParametersPointer);
}


// Get the parameters
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
const 
typename BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::ParametersType &
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::GetFixedParameters( void ) const
{
  RegionType resRegion = this->GetGridRegion(  );
  
  for (unsigned int i=0;i<NDimensions;i++)
    {
    this->m_FixedParameters[i] = (resRegion.GetSize())[i];
    }
  for (unsigned int i=0;i<NDimensions;i++)
    {
    this->m_FixedParameters[NDimensions+i] = (this->GetGridOrigin())[i];
    } 
  for (unsigned int i=0;i<NDimensions;i++)
    {
    this->m_FixedParameters[2*NDimensions+i] =  (this->GetGridSpacing())[i];
    }
  
  return (this->m_FixedParameters);
}



 
// Set the B-Spline coefficients using input images
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void 
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::SetCoefficientImage( ImagePointer images[] )
{
  if ( images[0] )
    {
    this->SetGridRegion( images[0]->GetBufferedRegion() );
    this->SetGridSpacing( images[0]->GetSpacing() );
    this->SetGridOrigin( images[0]->GetOrigin() );

    for( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_CoefficientImage[j] = images[j];
      }

    // Clean up buffered parameters
    m_InternalParametersBuffer = ParametersType( 0 );
    m_InputParametersPointer  = NULL;

    }

}  


// Print self
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::PrintSelf(std::ostream &os, Indent indent) const
{

  unsigned int j;

  this->Superclass::PrintSelf(os,indent);

  os << indent << "GridRegion: " << m_GridRegion << std::endl;
  os << indent << "GridOrigin: " << m_GridOrigin << std::endl;
  os << indent << "GridSpacing: " << m_GridSpacing << std::endl;

  os << indent << "CoefficientImage: [ ";
  for ( j = 0; j < SpaceDimension - 1; j++ )
    {
    os << m_CoefficientImage[j].GetPointer() << ", ";
    }
  os << m_CoefficientImage[j].GetPointer() << " ]" << std::endl;

  os << indent << "WrappedImage: [ ";
  for ( j = 0; j < SpaceDimension - 1; j++ )
    {
    os << m_WrappedImage[j].GetPointer() << ", ";
    }
  os << m_WrappedImage[j].GetPointer() << " ]" << std::endl;
 
  os << indent << "InputParametersPointer: " 
     << m_InputParametersPointer << std::endl;
  os << indent << "ValidRegion: " << m_ValidRegion << std::endl;
  os << indent << "LastJacobianIndex: " << m_LastJacobianIndex << std::endl;
  os << indent << "BulkTransform: ";
  os << m_BulkTransform.GetPointer() << std::endl;
  os << indent << "WeightsFunction: ";
  os << m_WeightsFunction.GetPointer() << std::endl;

  if ( m_BulkTransform )
    {
    os << indent << "BulkTransformType: " 
       << m_BulkTransform->GetNameOfClass() << std::endl;
    }
     
}

// Transform a point
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
bool 
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::InsideValidRegion( 
  const ContinuousIndexType& index ) const
{
  bool inside = true;

  if ( !m_ValidRegion.IsInside( index ) )
    {
    inside = false;
    }

  if ( inside && m_SplineOrderOdd )
    {
    typedef typename ContinuousIndexType::ValueType ValueType;
    for( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      if ( index[j] >= static_cast<ValueType>( m_ValidRegionLast[j] ) )
        { 
        inside = false;
        break;
        }
      }
    }

  return inside;
}



// Transform a point
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void 
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::TransformPoint( 
  const InputPointType & point, 
  OutputPointType & outputPoint, 
  WeightsType & weights ) const
{

  unsigned int j;
  IndexType supportIndex;

  InputPointType transformedPoint;
  if ( m_BulkTransform )
    {
    transformedPoint = m_BulkTransform->TransformPoint( point );
    }
  else
    {
    transformedPoint = point;
    }

  if ( m_CoefficientImage[0] )
    {

    ContinuousIndexType index;
    for ( j = 0; j < SpaceDimension; j++ )
      {
      index[j] = ( point[j] - m_GridOrigin[j] ) / m_GridSpacing[j];
      }

    // NOTE: if the support region does not lie totally within the grid
    // we assume zero displacement and return the input point
    bool inside = this->InsideValidRegion( index );
    if ( !inside )
      {
      outputPoint = transformedPoint;
      return;
      }

    // Compute interpolation weights
    m_WeightsFunction->Evaluate( index, weights, supportIndex );

    // For each dimension, correlate coefficient with weights
    RegionType supportRegion;
    supportRegion.SetSize( m_SupportSize );
    supportRegion.SetIndex( supportIndex );

    outputPoint.Fill( NumericTraits<ScalarType>::Zero );

    typedef ImageRegionConstIterator<ImageType> IteratorType;
    IteratorType m_Iterator[ SpaceDimension ];
    unsigned long counter = 0;

    for ( j = 0; j < SpaceDimension; j++ )
      {
      m_Iterator[j] = IteratorType( m_CoefficientImage[j], supportRegion );
      }

    while ( ! m_Iterator[0].IsAtEnd() )
      {

      // multiply weigth with coefficient
      for ( j = 0; j < SpaceDimension; j++ )
        {
        outputPoint[j] += static_cast<ScalarType>( 
          weights[counter] * m_Iterator[j].Get());
        }

      // go to next coefficient in the support region
      ++ counter;
      for ( j = 0; j < SpaceDimension; j++ )
        {
        ++( m_Iterator[j] );
        }
      }
  
    // return results
    for ( j = 0; j < SpaceDimension; j++ )
      {
      outputPoint[j] += transformedPoint[j];
      }

    }
    else
    {

    itkWarningMacro( << "B-spline coefficients have not been set" );

    for ( j = 0; j < SpaceDimension; j++ )
      {
      outputPoint[j] = transformedPoint[j];
      }

    }

}


// Transform a point
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
typename BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::OutputPointType
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::TransformPoint(const InputPointType &point) const 
{
  
  WeightsType weights( m_WeightsFunction->GetNumberOfWeights() );
  OutputPointType outputPoint;

  this->TransformPoint( point, outputPoint, weights );

  return outputPoint;

}

 
// Compute the Jacobian in one position 
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
const 
typename BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::JacobianType & 
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::GetJacobian( const InputPointType & point ) const
{
  // Can only compute Jacobian if parameters are set via
  // SetParameters or SetParametersByValue
  if( m_InputParametersPointer == NULL )
    {
    itkExceptionMacro( <<"Cannot compute Jacobian: parameters not set" );
    }

  // Zero all components of jacobian
  // NOTE: for efficiency, we only need to zero out the coefficients
  // that got fill last time this method was called.
  RegionType supportRegion;
  supportRegion.SetSize( m_SupportSize );
  supportRegion.SetIndex( m_LastJacobianIndex );

  typedef ImageRegionIterator<JacobianImageType> IteratorType;
  IteratorType m_Iterator[ SpaceDimension ];
  unsigned int j;

  for ( j = 0; j < SpaceDimension; j++ )
    {
    m_Iterator[j] = IteratorType( m_JacobianImage[j], supportRegion );
    }

  while ( ! m_Iterator[0].IsAtEnd() )
    {

    // zero out jacobian elements
    for ( j = 0; j < SpaceDimension; j++ )
      {
      m_Iterator[j].Set( NumericTraits<JacobianPixelType>::Zero );
      }

    for ( j = 0; j < SpaceDimension; j++ )
      {
      ++( m_Iterator[j] );
      }
    }

 
  ContinuousIndexType index;
  for ( j = 0; j < SpaceDimension; j++ )
    {
    index[j] = ( point[j] - m_GridOrigin[j] ) / m_GridSpacing[j];
    }

  // NOTE: if the support region does not lie totally within the grid
  // we assume zero displacement and return the input point
  if ( !this->InsideValidRegion( index ) )
    {
    return this->m_Jacobian;
    }

  // Compute interpolation weights
  WeightsType weights( m_WeightsFunction->GetNumberOfWeights() );
  IndexType supportIndex;

  m_WeightsFunction->Evaluate( index, weights, supportIndex );
  m_LastJacobianIndex = supportIndex;

  // For each dimension, copy the weight to the support region
  supportRegion.SetIndex( supportIndex );
  unsigned long counter = 0;

  for ( j = 0; j < SpaceDimension; j++ )
    {
    m_Iterator[j] = IteratorType( m_JacobianImage[j], supportRegion );
    }

  while ( ! m_Iterator[0].IsAtEnd() )
    {

    // copy weight to jacobian image
    for ( j = 0; j < SpaceDimension; j++ )
      {
      m_Iterator[j].Set( static_cast<JacobianPixelType>( weights[counter] ) );
      }

    // go to next coefficient in the support region
    ++ counter;
    for ( j = 0; j < SpaceDimension; j++ )
      {
      ++( m_Iterator[j] );
      }
    }


  // Return the results
  return this->m_Jacobian;
}



// Compute the Jacobian in one position 
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::GetJacobian( const InputPointType & point, ParameterIndexArrayType& indexes, WeightsType& weights )
{

  RegionType supportRegion;
  supportRegion.SetSize( m_SupportSize );
  const PixelType * basePointer = m_CoefficientImage[0]->GetBufferPointer();

  unsigned int j;

  ContinuousIndexType index;
  for ( j = 0; j < SpaceDimension; j++ )
  {
    index[j] = ( point[j] - m_GridOrigin[j] ) / m_GridSpacing[j];
  }

  // NOTE: if the support region does not lie totally within the grid
  // we assume zero displacement and return the input point
  if ( !this->InsideValidRegion( index ) )
  {
    weights.Fill(0.0);
    indexes.Fill(0);
    return;
  }
  
  // Compute interpolation weights
  IndexType supportIndex;

  m_WeightsFunction->Evaluate( index, weights, supportIndex );

  // For each dimension, copy the weight to the support region
  supportRegion.SetIndex( supportIndex );
  unsigned long counter = 0;

  typedef ImageRegionIterator<JacobianImageType> IteratorType;

  IteratorType m_Iterator = IteratorType( m_CoefficientImage[0], supportRegion );


  while ( ! m_Iterator.IsAtEnd() )
    {


    indexes[counter] = &(m_Iterator.Value()) - basePointer;

    // go to next coefficient in the support region
    ++ counter;
    ++m_Iterator;
    
    }

}

template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
unsigned int 
BSplineDeformableTransformOpt<TScalarType, NDimensions,VSplineOrder>
::GetNumberOfAffectedWeights() const
{
  return m_WeightsFunction->GetNumberOfWeights();
}


} // namespace

#endif
