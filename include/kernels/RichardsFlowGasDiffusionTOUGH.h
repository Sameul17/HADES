//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/**
 * Convective flux of component k in fluid phase alpha.
 * A fully-updwinded version is implemented, where the mobility
 * of the upstream nodes is used.
 */
class RichardsFlowGasDiffusionTOUGH : public Kernel
{
public:
  static InputParameters validParams();

  RichardsFlowGasDiffusionTOUGH(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Temperature of each quadrature point
  const VariableValue & _T;
  const VariableGradient & _grad_T;
  unsigned _T_id;

  // Water density of each quadrature point
  const MaterialProperty<Real> & _Wdensity;

  // Bentonite porosity of each quadrature point
  const MaterialProperty<Real> & _porosity;

  // Van genuchten parameter of each quadrature point
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _lambda;

  /// Gas pressure of each quadrature point
//  const VariableValue & _Pgas;

};
