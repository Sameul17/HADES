//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowAdvectiveFlux.h"

registerMooseObject("HadesApp", RichardsFlowAdvectiveFlux);

InputParameters
RichardsFlowAdvectiveFlux::validParams()
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

RichardsFlowAdvectiveFlux::RichardsFlowAdvectiveFlux(const InputParameters & parameters)
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
RichardsFlowAdvectiveFlux::computeQpResidual()
{
  Real rho_w = 1000;

  return _grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_u[_qp] - rho_w * _gravity);
}

Real
RichardsFlowAdvectiveFlux::computeQpJacobian()
{

  Real rho_w = 1000;

  return _grad_test[_i][_qp] * rho_w * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * _grad_phi[_j][_qp];
}

Real
RichardsFlowAdvectiveFlux::computeQpOffDiagJacobian(unsigned int jvar)
{

  return 0.0;
}
