//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowMassTimeDerivative.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include <limits>

registerMooseObject("HadesApp", RichardsFlowMassTimeDerivative);

InputParameters
RichardsFlowMassTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point"); // Temperature [K]
  params.addCoupledVar("Pgas", 1E5, "The gas pressure of each quadrature point"); // Gas pressure [Pa]
  params.addRequiredParam<MaterialPropertyName>("porosity", "The porostiy of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The density of water of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The parameter of Van Genuchten equation of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The parameter of Van Genuchten equation of each quadrature point");

  params.addClassDescription("Derivative of fluid-component mass with respect to time. The variable u is capillary pressure in here");
  return params;
}

RichardsFlowMassTimeDerivative::RichardsFlowMassTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
   _u_old(valueOld()),
   _T(coupledValue("Temp")),
   _T_old(coupledValue("Temp")),
   _T_id(coupled("Temp")),
   _Pgas(coupledValue("Pgas")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _porosity(getMaterialProperty<Real>("porosity"))
{
}

Real
RichardsFlowMassTimeDerivative::computeQpResidual()
{
  Real Rv = 461.5; // Specific gas constant [J/(kgK)] Gas constant / Molar mass of water
  Real rho_w = 1000; // Water density

  // The gas density given by the Clapeyron equation and Dalton's law
  Real gas_density = 0.001 * exp(19.84 - 4976 / _T[_qp] - _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real gas_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] - _u_old[_qp] / (rho_w * Rv * _T_old[_qp]));
 
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _u_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Time derivative terms
  Real dSat_dt = (Sat - Sat_old) / _dt;
  Real dgas_density_dt = (gas_density - gas_density_old) / _dt;

  return _test[_i][_qp] * _porosity[_qp] * ((rho_w - gas_density) * dSat_dt + (1 - Sat) * dgas_density_dt);
}

Real
RichardsFlowMassTimeDerivative::computeQpJacobian()
{
  Real Rv = 461.5; // Specific gas constant [J/(kgK)] Gas constant / Molar mass of water
  Real rho_w = 1000; // Water density

  // The gas density given by the Clapeyron equation and Dalton's law
  Real gas_density = 0.001 * exp(19.84 - 4976 / _T[_qp] - _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real gas_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] - _u_old[_qp] / (rho_w * Rv * _T_old[_qp]));

  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _u_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Time derivative terms
  Real dT_dt = (_T[_qp] - _T_old[_qp]) / _dt;
  Real du_dt = (_u[_qp] - _u_old[_qp]) / _dt;
  Real dSat_dt = (Sat - Sat_old) / _dt;
  Real dgas_density_dt = (gas_density - gas_density_old) / _dt;

  // Jacobian term 
  Real dgas_density_dPcap = -gas_density / (rho_w * Rv * _T[_qp]) * _phi[_j][_qp];
  Real dSat_dPcap = _lambda[_qp] * _alpha[_qp] / (1 - _lambda[_qp]) * pow((-_alpha[_qp] * _u[_qp]),(_lambda[_qp] / (1 - _lambda[_qp]))) * pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1/(1 - _lambda[_qp])))),(-_lambda[_qp] - 1)) * _phi[_j][_qp];

  Real dm = -_test[_i][_qp] * _porosity[_qp] * dSat_dt * dgas_density_dPcap;
  dm += _test[_i][_qp] * _porosity[_qp] * (rho_w - gas_density) * dSat_dPcap / _dt;
  dm -= _test[_i][_qp] * _porosity[_qp] * dSat_dPcap * dgas_density_dt;
  dm += _test[_i][_qp] * _porosity[_qp] * (1 - Sat) * dgas_density_dPcap / _dt;

  return dm;
}

Real
RichardsFlowMassTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real Rv = 461.5; // Specific gas constant [J/(kgK)] Gas constant / Molar mass of water
  Real rho_w = 1000; // Water density

  // The gas density given by the Clapeyron equation and Dalton's law
  Real gas_density = 0.001 * exp(19.84 - 4976 / _T[_qp] - _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real gas_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] - _u_old[_qp] / (rho_w * Rv * _T_old[_qp]));

  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _u_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Time derivative terms
  Real dSat_dt = (Sat - Sat_old) / _dt;
  Real dT_dt = (_T[_qp] - _T_old[_qp]) / _dt;
  Real du_dt = (_u[_qp] - _u_old[_qp]) / _dt;
  Real dgas_density_dt = (gas_density - gas_density_old) / _dt;

  // Jacobian term
  Real dphi_dt = _phi[_j][_qp] / _dt;
  Real dgas_density_dT = gas_density * (4976 / pow(_T[_qp],2) + _u[_qp] / (rho_w * Rv * pow(_T[_qp],2))) * _phi[_j][_qp];

  Real dm = -_test[_i][_qp] * _porosity[_qp] * dgas_density_dT * dSat_dt;
      dm += _test[_i][_qp] * _porosity[_qp] * (1 - Sat) * dgas_density_dT / _dt;

  if (jvar == _T_id)
	return dm;
  else
	return 0.0;  
}

