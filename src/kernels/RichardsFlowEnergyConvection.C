//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowEnergyConvection.h"

registerMooseObject("HadesApp", RichardsFlowEnergyConvection);

InputParameters
RichardsFlowEnergyConvection::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredCoupledVar("Pwater", "The temperature of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Intrinsic_permeability", "The intrinsic permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Relative_permeability", "The Relative permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_viscosity", "The water viscosity of each quadrature point");
  params.addRequiredParam<RealVectorValue>("gravity", "Gravitational acceleration vector downwards (m/s^2)");

  params.addClassDescription("Calculate energy convection in porous medium");
  return params;
}

RichardsFlowEnergyConvection::RichardsFlowEnergyConvection(const InputParameters & parameters)
  : Kernel(parameters),
   _Pwater(coupledValue("Pwater")),
   _grad_Pwater(coupledGradient("Pwater")),
   _Pwater_id(coupled("Pwater")),
   _IPerm(getMaterialProperty<Real>("Intrinsic_permeability")),
   _RPerm(getMaterialProperty<Real>("Relative_permeability")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _Wviscosity(getMaterialProperty<Real>("Water_viscosity")),
   _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
RichardsFlowEnergyConvection::computeQpResidual()
{
  Real rho_w = 1E3;

  return -_grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_Pwater[_qp] - rho_w * _gravity) * _u[_qp];
//  return -_grad_test[_i][_qp] * _Wdensity[_qp] * _Cp[_qp] * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_Pwater[_qp] - _Wdensity[_qp] * _gravity) * _u[_qp];
}

Real
RichardsFlowEnergyConvection::computeQpJacobian()
{
  Real rho_w = 1E3;
 
//  return 0.0;
  return -_grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_Pwater[_qp] - rho_w * _gravity) * _phi[_j][_qp];
//  return -_grad_test[_i][_qp] * _Wdensity[_qp] * _Cp[_qp] * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_Pwater[_qp] - _Wdensity[_qp] * _gravity) * _phi[_j][_qp];
}

Real
RichardsFlowEnergyConvection::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real rho_w = 1E3;

  if (jvar == _Pwater_id)
    return -_grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_phi[_j][_qp]) * _u[_qp];
//	  return -_grad_test[_i][_qp] * _Wdensity[_qp] * _Cp[_qp] * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_phi[_j][_qp]) * _u[_qp];
  else
	  return 0.0;
}
