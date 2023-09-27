//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowEnergyTimeDerivative.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include <limits>

registerMooseObject("HadesApp", RichardsFlowEnergyTimeDerivative);

InputParameters
RichardsFlowEnergyTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addCoupledVar("Pwater", 1E5, "The capillary pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Medium_density", "The density of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Medium_specific_heat", "The specific heat of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The density of water of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Porosity", "The porosity of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The lambda of Van Genuchten equation of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The alpha of Van Genuchten equation of each quadrature point");
  

  params.addClassDescription("Derivative of temperature with respect to time. The variable u is temperature in here");
  return params;
}

RichardsFlowEnergyTimeDerivative::RichardsFlowEnergyTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
   _u_old(valueOld()),
   _P(coupledValue("Pwater")),
   _MDensity(getMaterialProperty<Real>("Medium_density")),
   _MCp(getMaterialProperty<Real>("Medium_specific_heat")),
   _WDensity(getMaterialProperty<Real>("Water_density")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _porosity(getMaterialProperty<Real>("Porosity"))
{
}

Real
RichardsFlowEnergyTimeDerivative::computeQpResidual()
{
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real WCp = 4.2E3; // Water specific heat [J/kgk]
 
  return _test[_i][_qp] * ((1 - _porosity[_qp]) * _MDensity[_qp] * _MCp[_qp] + Sat * _porosity[_qp] * _WDensity[_qp] * WCp) * (_u[_qp] - _u_old[_qp]) / _dt;
}

Real
RichardsFlowEnergyTimeDerivative::computeQpJacobian()
{
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real WCp = 4.2E3; // Water specific heat [J/kgk]

  return _test[_i][_qp] * ((1 - _porosity[_qp]) * _MDensity[_qp] * _MCp[_qp] + Sat * _porosity[_qp] * _WDensity[_qp] * WCp) * _phi[_j][_qp] / _dt;
}

Real
RichardsFlowEnergyTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real ddSat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp] - 1));
  Real dSat_dVar = _lambda[_qp] * _alpha[_qp] / (1 - _lambda[_qp]) * pow((-_alpha[_qp] * _P[_qp]),(_lambda[_qp] / (1 - _lambda[_qp]))) * ddSat * _phi[_j][_qp];
  Real WCp = 4.2E3; // Water specific heat [J/kgk]

  return _test[_i][_qp] * dSat_dVar * _WDensity[_qp] * WCp * (_u[_qp] - _u_old[_qp]) / _dt;
}

