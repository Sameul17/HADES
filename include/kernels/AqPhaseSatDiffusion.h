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


class AqPhaseSatDiffusion : public Kernel
{
public:
  static InputParameters validParams();

  AqPhaseSatDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  ///Coupled Variables
  const VariableValue & _T;
  unsigned _T_id;

  const VariableValue & _P;
  unsigned _P_id;

  /// Porosity at the nodes, but it can depend on grad(variables) which are actually evaluated at the qps
  const MaterialProperty<Real> & _porosity;

  /// Tortuosity at the nodess
  const MaterialProperty<Real> & _tortuosity;

  /// Calculation of diffusion coefficient
  const MaterialProperty<Real> & _Diffusivity;
//  const MaterialProperty<Real> & _D_aq;
//  const MaterialProperty<Real> & _EA_aq;

  /// Van Genuchten parameter
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _lambda;

};
