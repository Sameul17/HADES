//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowEnergyTimeDerivativeSol.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include <limits>

registerMooseObject("HadesApp", RichardsFlowEnergyTimeDerivativeSol);

InputParameters
RichardsFlowEnergyTimeDerivativeSol::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("Medium_density", "The density of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Medium_specific_heat", "The specific heat of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The density of water of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Porosity", "The porosity of medium of each quadrature point");
  

  params.addClassDescription("Derivative of temperature with respect to time. The variable u is temperature in here");
  return params;
}

RichardsFlowEnergyTimeDerivativeSol::RichardsFlowEnergyTimeDerivativeSol(const InputParameters & parameters)
  : TimeKernel(parameters),
   _u_old(valueOld()),
   _MDensity(getMaterialProperty<Real>("Medium_density")),
   _MCp(getMaterialProperty<Real>("Medium_specific_heat")),
   _WDensity(getMaterialProperty<Real>("Water_density")),
   _porosity(getMaterialProperty<Real>("Porosity"))
{
}

Real
RichardsFlowEnergyTimeDerivativeSol::computeQpResidual()
{
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real WCp = 4.2E3; // Water specific heat [J/kgk]
 
  return _test[_i][_qp] * (_MDensity[_qp] * _MCp[_qp]) * (_u[_qp] - _u_old[_qp]) / _dt;
}

Real
RichardsFlowEnergyTimeDerivativeSol::computeQpJacobian()
{
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  
  Real WCp = 4.2E3; // Water specific heat [J/kgk]

  return _test[_i][_qp] * (_MDensity[_qp] * _MCp[_qp]) * _phi[_j][_qp] / _dt;
}

Real
RichardsFlowEnergyTimeDerivativeSol::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real WCp = 4.2E3; // Water specific heat [J/kgk]

  return 0;
}

