/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkNormalizedCorrelationTwoImageToOneImageMetric.txx,v $
Language:  C++
Date:      $Date: 2010/12/20 $
Version:   $Revision: 1.0 $
Author:    Jian Wu (eewujian@hotmail.com)
Univerisity of Florida
Virginia Commonwealth University

This program was modified from the ITK program--itkNormalizedCorrelationImageToImageMetric.txx

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkNormalizedCorrelationTwoImageToOneImageMetric_txx
#define _itkNormalizedCorrelationTwoImageToOneImageMetric_txx

#include "itkNormalizedCorrelationTwoImageToOneImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

	/*
	 * Constructor
	 */
	template <class TFixedImage, class TMovingImage> 
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,TMovingImage>
		::NormalizedCorrelationTwoImageToOneImageMetric()
		{
			m_SubtractMean = false;
		}

	/**
	 * Initialize
	 */

	template <class TFixedImage, class TMovingImage>
		void NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage, TMovingImage>
		::Initialize(void) throw (ExceptionObject)
		{
			//Initialize the base class
			Superclass::Initialize();


		}


	/*
	 * Get the match Measure
	 */
	template <class TFixedImage, class TMovingImage> 
		typename NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,TMovingImage>::MeasureType
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,TMovingImage>
		::GetValue( const TransformParametersType & parameters ) const
		{

			FixedImageConstPointer fixedImage = this->m_FixedImage;

			if( !fixedImage ) 
			{
				itkExceptionMacro( << "Fixed image1 has not been assigned" );
			}


			typedef  itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedIteratorType;

			typedef  typename NumericTraits< MeasureType >::AccumulateType AccumulateType;


			// Calculate the measure value between fixed image 1 and the moving image

			FixedIteratorType ti1( fixedImage, this->GetFixedImageRegion() );

			typename FixedImageType::IndexType index;

			MeasureType measure1;

			this->m_NumberOfPixelsCounted = 0;

			this->SetTransformParameters( parameters );

			AccumulateType sff = NumericTraits< AccumulateType >::Zero;
			AccumulateType smm = NumericTraits< AccumulateType >::Zero;
			AccumulateType sfm = NumericTraits< AccumulateType >::Zero;
			AccumulateType sf  = NumericTraits< AccumulateType >::Zero;
			AccumulateType sm  = NumericTraits< AccumulateType >::Zero;
			typename Superclass::InputPointType inputPoint;

			while(!ti1.IsAtEnd())
			{

				index = ti1.GetIndex();

				fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

				if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
				{
					++ti1;
					continue;
				}

				//    typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

				if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( inputPoint ) )
				{
					++ti1;
					continue;
				}

				if( this->m_Interpolator->IsInsideBuffer( inputPoint ) )
				{
					//uso del 2er interpolar evaluando solo un punto actual de toda la imagen fija
					const RealType movingValue  = this->m_Interpolator->Evaluate( inputPoint );
					const RealType fixedValue   = ti1.Get();
					sff += fixedValue  * fixedValue;
					smm += movingValue * movingValue;
					sfm += fixedValue  * movingValue;
					if ( this->m_SubtractMean )
					{
						sf += fixedValue;
						sm += movingValue;
					}
					this->m_NumberOfPixelsCounted++;
				}

				++ti1;
			}

			if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
			{
				sff -= ( sf * sf / this->m_NumberOfPixelsCounted );
				smm -= ( sm * sm / this->m_NumberOfPixelsCounted );
				sfm -= ( sf * sm / this->m_NumberOfPixelsCounted );
			}

			RealType denom = -1.0 * sqrt( sff * smm );

			if( this->m_NumberOfPixelsCounted > 0 && denom != 0.0)
			{
				measure1 = sfm / denom;
			}
			else
			{
				measure1 = NumericTraits< MeasureType >::Zero;
			}



			return measure1;

		}





	/*
	 * Get the Derivative Measure
	 */
	template < class TFixedImage, class TMovingImage> 
		void
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,TMovingImage>
		::GetDerivative( const TransformParametersType & parameters,
				DerivativeType & derivative ) const
		{
			// under construction
			//no contiene derivada basada en diferencias finitas
			//asi que la metrica es solo retornado como valor numerico.
		}


	/*
	 * Get both the match Measure and theDerivative Measure 
	 */
	template <class TFixedImage, class TMovingImage> 
		void
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,TMovingImage>
		::GetValueAndDerivative(const TransformParametersType & parameters, 
				MeasureType & value, DerivativeType  & derivative) const
		{
			// under construction
		}

	template < class TFixedImage, class TMovingImage> 
		void
		NormalizedCorrelationTwoImageToOneImageMetric<TFixedImage,TMovingImage>
		::PrintSelf(std::ostream& os, Indent indent) const
		{
			Superclass::PrintSelf(os, indent);
			os << indent << "SubtractMean: " << m_SubtractMean << std::endl;
		}

} // end namespace itk


#endif
