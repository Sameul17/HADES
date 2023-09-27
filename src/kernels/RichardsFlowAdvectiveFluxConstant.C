//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowAdvectiveFluxConstant.h"

registerMooseObject("HadesApp", RichardsFlowAdvectiveFluxConstant);

InputParameters
RichardsFlowAdvectiveFluxConstant::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point");
//  params.addRequiredCoupledVar("Pgas", "The gas pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Intrinsic_permeability", "The intrinsic permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Relative_permeability", "The Relative permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_viscosity", "The dynamic viscosity of water of each quadrature point");
  params.addRequiredParam<RealVectorValue>("gravity","Gravitational acceleration vector downwards (m/s^2)");

  params.addClassDescription("Calculate water advection in porous medium");
  return params;
}

RichardsFlowAdvectiveFluxConstant::RichardsFlowAdvectiveFluxConstant(const InputParameters & parameters)
  : Kernel(parameters),
   _T(coupledValue("Temp")),
//   _Pgas(coupledValue("Pgas")),
   _IPerm(getMaterialProperty<Real>("Intrinsic_permeability")),
   _RPerm(getMaterialProperty<Real>("Relative_permeability")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _Wviscosity(getMaterialProperty<Real>("Water_viscosity")),
   _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
RichardsFlowAdvectiveFluxConstant::computeQpResidual()
{

  return _grad_test[_i][_qp] * _Wdensity[_qp] * _IPerm[_qp] * _RPerm[_qp] / 3.17652E-11 * (_grad_u[_qp] - _Wdensity[_qp] * _gravity);
}

Real
RichardsFlowAdvectiveFluxConstant::computeQpJacobian()
{

  return _grad_test[_i][_qp] * _Wdensity[_qp] * _IPerm[_qp] * _RPerm[_qp] / 3.17652E-11 * _grad_phi[_j][_qp];
}

Real
RichardsFlowAdvectiveFluxConstant::computeQpOffDiagJacobian(unsigned int jvar)
{

  return 0.0;
}
