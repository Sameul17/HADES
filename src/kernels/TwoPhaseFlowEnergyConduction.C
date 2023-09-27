//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowEnergyConduction.h"

registerMooseObject("HadesApp", TwoPhaseFlowEnergyConduction);

InputParameters
TwoPhaseFlowEnergyConduction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("Thermal_conductivity", "The effective thermal conductivity of each quadrature point");

  params.addClassDescription("Calculate the Gas Diffusion behavior in porous medium");
  return params;
}

TwoPhaseFlowEnergyConduction::TwoPhaseFlowEnergyConduction(const InputParameters & parameters)
  : Kernel(parameters),
  _thermal_conductivity(getMaterialProperty<Real>("Thermal_conductivity"))
{
}

Real
TwoPhaseFlowEnergyConduction::computeQpResidual()
{
  return _grad_test[_i][_qp] * _thermal_conductivity[_qp] * _grad_u[_qp];
}

Real
TwoPhaseFlowEnergyConduction::computeQpJacobian()
{
  return _grad_test[_i][_qp] * _thermal_conductivity[_qp] * _grad_phi[_j][_qp];
}

Real
TwoPhaseFlowEnergyConduction::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0;
}
