//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowEnergyConduction.h"

registerMooseObject("HadesApp", RichardsFlowEnergyConduction);

InputParameters
RichardsFlowEnergyConduction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("Thermal_conductivity", "The effective thermal conductivity of each quadrature point");

  params.addClassDescription("Calculate the Gas Diffusion behavior in porous medium");
  return params;
}

RichardsFlowEnergyConduction::RichardsFlowEnergyConduction(const InputParameters & parameters)
  : Kernel(parameters),
  _thermal_conductivity(getMaterialProperty<Real>("Thermal_conductivity"))
{
}

Real
RichardsFlowEnergyConduction::computeQpResidual()
{
  return _grad_test[_i][_qp] * _thermal_conductivity[_qp] * _grad_u[_qp];
}

Real
RichardsFlowEnergyConduction::computeQpJacobian()
{
  return _grad_test[_i][_qp] * _thermal_conductivity[_qp] * _grad_phi[_j][_qp];
}

Real
RichardsFlowEnergyConduction::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
