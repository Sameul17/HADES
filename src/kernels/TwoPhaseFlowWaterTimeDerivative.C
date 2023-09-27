//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowWaterTimeDerivative.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include <limits>

registerMooseObject("HadesApp", TwoPhaseFlowWaterTimeDerivative);

InputParameters
TwoPhaseFlowWaterTimeDerivative::validParams()
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

TwoPhaseFlowWaterTimeDerivative::TwoPhaseFlowWaterTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
   _u_old(valueOld()),
   _T(coupledValue("Temp")),
   _T_old(coupledValue("Temp")),
   _T_id(coupled("Temp")),
   _Pgas(coupledValue("Pgas")),
   _Pgas_old(coupledValue("Pgas")),
   _Pgas_id(coupled("Pgas")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _porosity(getMaterialProperty<Real>("porosity"))
{
}

  // Primary variable: Capillary pressure (Pcap)

Real
TwoPhaseFlowWaterTimeDerivative::computeQpResidual()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real rho_w = 1000; // Water density

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] + _u_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _u_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Time derivative termsi
  Real dSat_dt = (Sat - Sat_old) / _dt;
  Real dvapor_density_dt = (vapor_density - vapor_density_old) / _dt;

  return _test[_i][_qp] * _porosity[_qp] * ((rho_w - vapor_density) * dSat_dt + (1 - Sat) * dvapor_density_dt);
}

Real
TwoPhaseFlowWaterTimeDerivative::computeQpJacobian()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real rho_w = 1000; // Water density

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] + _u_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _u_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Time derivative termsi
  Real dSat_dt = (Sat - Sat_old) / _dt;
  Real dvapor_density_dt = (vapor_density - vapor_density_old) / _dt;

  // Jacobian term, The variable is capillary presure
  Real dvapor_density_dPcap = Mw / (rho_w * R * _T[_qp]) * vapor_density * _phi[_j][_qp];
  Real dSat_dPcap = _alpha[_qp] * _lambda[_qp] / (1 - _lambda[_qp]) * pow((-_alpha[_qp] * _u[_qp]),(_lambda[_qp] / (1 - _lambda[_qp]))) * pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp] - 1)) * _phi[_j][_qp];

  // Jacobian term
  Real dm = -_test[_i][_qp] * _porosity[_qp] * dvapor_density_dPcap * dSat_dt;
      dm += _test[_i][_qp] * _porosity[_qp] * (rho_w - vapor_density) * dSat_dPcap / _dt;
      dm -= _test[_i][_qp] * _porosity[_qp] * dSat_dPcap * dvapor_density_dt;
      dm += _test[_i][_qp] * _porosity[_qp] * (1 - Sat) * dvapor_density_dPcap / _dt;
	  
  return dm;
//  return 0;
}

Real
TwoPhaseFlowWaterTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real rho_w = 1000; // Water density

  // The gas density given by the Clapeyron equation and Dalton's law
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] + _u_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _u_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Time derivative terms
  Real dSat_dt = (Sat - Sat_old) / _dt;
  Real dvapor_density_dt = (vapor_density - vapor_density_old) / _dt;

  // OffDiagJacobian term
  Real dvapor_density_dT = ((4976 / pow(_T[_qp],2)) - Mw * _u[_qp] / (rho_w * R * pow(_T[_qp],2))) * vapor_density * _phi[_j][_qp];

  // Jacobian term
  Real dm_T = -_test[_i][_qp] * _porosity[_qp] * dvapor_density_dT * dSat_dt;
      dm_T += _test[_i][_qp] * _porosity[_qp] * (1 - Sat) * dvapor_density_dT / _dt;
  Real dm_Pgas = 0;

  if (jvar == _T_id)
	return dm_T;
  else if (jvar == _Pgas_id)
	return dm_Pgas;
  else
        return 0.0;  
}

