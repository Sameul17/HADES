//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowLiquidAdvectiveFlux.h"

registerMooseObject("HadesApp", TwoPhaseFlowLiquidAdvectiveFlux);

InputParameters
TwoPhaseFlowLiquidAdvectiveFlux::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Pgas", 1E5, "The temperature of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Intrinsic_permeability", "The intrinsic permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Relative_permeability", "The Relative permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_viscosity", "The dynamic viscosity of water of each quadrature point");
  params.addRequiredParam<RealVectorValue>("gravity","Gravitational acceleration vector downwards (m/s^2)");

  params.addClassDescription("Calculate water advection in porous medium");
  return params;
}

TwoPhaseFlowLiquidAdvectiveFlux::TwoPhaseFlowLiquidAdvectiveFlux(const InputParameters & parameters)
  : Kernel(parameters),
   _Pgas(coupledValue("Pgas")),
   _Pgas_old(coupledValue("Pgas")),
   _grad_Pgas(coupledGradient("Pgas")),
   _Pgas_id(coupled("Pgas")),
   _IPerm(getMaterialProperty<Real>("Intrinsic_permeability")),
   _RPerm(getMaterialProperty<Real>("Relative_permeability")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _Wviscosity(getMaterialProperty<Real>("Water_viscosity")),
   _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
TwoPhaseFlowLiquidAdvectiveFlux::computeQpResidual()
{
  Real rho_w = 1E3;

  return _grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (-_grad_Pgas[_qp] + _grad_u[_qp] - rho_w * _gravity);
//  return _grad_test[_i][_qp] * _Wdensity[_qp] * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_Pgas[_qp] - _grad_u[_qp] - _Wdensity[_qp] * _gravity);
}

Real
TwoPhaseFlowLiquidAdvectiveFlux::computeQpJacobian()
{
  Real rho_w = 1E3;

  return _grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * _grad_phi[_j][_qp];
//  return -_grad_test[_i][_qp] * _Wdensity[_qp] * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * _grad_phi[_j][_qp];
}

Real
TwoPhaseFlowLiquidAdvectiveFlux::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real rho_w = 1E3;

  if (jvar == _Pgas_id)
       return -_grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * _grad_phi[_j][_qp];
//       return _grad_test[_i][_qp] * _Wdensity[_qp] * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * _grad_phi[_j][_qp];
  else
       return 0.0;
}
