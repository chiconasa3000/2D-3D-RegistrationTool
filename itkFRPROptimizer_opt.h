/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkFRPROptimizer_opt_h
#define itkFRPROptimizer_opt_h

#include "itkPowellOptimizer.h"
#include "ITKOptimizersExport.h"

namespace itk
{
/** \class FRPROptimizer
 * \brief Implements Fletch-Reeves & Polak-Ribiere optimization using dBrent
 * line search.
 *
 * This optimizer needs a cost function.
 * This optimizer needs to be able to compute partial derivatives of the
 *    cost function with respect to each parameter.
 *
 * The SetStepLength determines the initial distance to step in a line
 * direction when bounding the minimum (using bracketing triple spaced
 * using a derivative-based search strategy).
 *
 * The StepTolerance terminates optimization when the parameter values are
 * known to be within this (scaled) distance of the local extreme.
 *
 * The ValueTolerance terminates optimization when the cost function values at
 * the current parameters and at the local extreme are likely (within a second
 * order approximation) to be within this is tolerance.
 *
 * \ingroup Numerics Optimizers
 *
 * \ingroup ITKOptimizers
 */
//// CUDA Runtime
//#include <cuda_runtime.h>
//
//// Using updated (v2) interfaces for CUBLAS and CUSPARSE
//#include <cusparse.h>
//#include <cublas_v2.h>
//
//// Utilities and system includes
//#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
//#include <helper_cuda.h>       // helper for CUDA error checking
//


class ITKOptimizers_EXPORT FRPROptimizer_opt:
  public PowellOptimizer
{
public:
  /** Standard "Self" type alias. */
  using Self = FRPROptimizer_opt;
  using Superclass = PowellOptimizer;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using ParametersType = SingleValuedNonLinearOptimizer::ParametersType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FRPROptimizer_opt, PowellOptimizer);

  /** Type of the Cost Function   */
  using CostFunctionType = SingleValuedCostFunction;
  using CostFunctionPointer = CostFunctionType::Pointer;

  /** Convert gradient to a unit length vector */
  itkSetMacro(UseUnitLengthGradient, bool);
  itkGetConstMacro(UseUnitLengthGradient, bool);

  /** Start optimization. */
  void StartOptimization() override;

  /** Set it to the Fletch-Reeves optimizer */
  void SetToFletchReeves();

  /** Set it to the Fletch-Reeves optimizer */
  void SetToPolakRibiere();

protected:
  FRPROptimizer_opt();
  ~FRPROptimizer_opt() override;

  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Get the value of the n-dimensional cost function at this scalar step
   * distance along the current line direction from the current line origin.
   * Line origin and distances are set via SetLine */
  virtual void GetValueAndDerivative(ParametersType & p, double *val,
                                     ParametersType *xi);

  virtual void   LineOptimize(ParametersType *p, ParametersType & xi,
                              double *val);

  virtual void   LineOptimize(ParametersType *p, ParametersType & xi,
                              double *val,
                              ParametersType & tempCoord);

private:
  FRPROptimizer_opt(const FRPROptimizer_opt &) = delete;

  typedef enum {
    FletchReeves,
    PolakRibiere
    }               OptimizationType;

  OptimizationType m_OptimizationType;

  bool m_UseUnitLengthGradient;
}; // end of class
} // end of namespace itk

#endif
