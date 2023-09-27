//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * Auxiliary kernel responsible for computing the components of the flux vector
 * in diffusion problems
 */
class PorousDiffusionFlux : public AuxKernel
{
public:
  static InputParameters validParams();

  PorousDiffusionFlux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  const VariableValue & _T;
  const VariableValue & _P;

  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _lambda;

  const MaterialProperty<Real> & _porosity;
  //
  /// Tortuosity at the nodess
  const MaterialProperty<Real> & _tortuosity;
  //
  //                    /// Calculation of diffusion coefficient
  const MaterialProperty<Real> & _D_aq;
  const MaterialProperty<Real> & _EA_aq;

  const MooseArray<Point> & _normals;
};
