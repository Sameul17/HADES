//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowAirTimeDerivative.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include <limits>

registerMooseObject("HadesApp", TwoPhaseFlowAirTimeDerivative);

InputParameters
TwoPhaseFlowAirTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point"); // Temperature [K]
  params.addCoupledVar("Pcap", 1E5, "The capillary pressure of each quadrature point"); // Gas pressure [Pa]
  params.addRequiredParam<MaterialPropertyName>("porosity", "The porostiy of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The density of water of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The parameter of Van Genuchten equation of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The parameter of Van Genuchten equation of each quadrature point");

  params.addClassDescription("Derivative of fluid-component mass with respect to time. The variable u is capillary pressure in here");
  return params;
}

TwoPhaseFlowAirTimeDerivative::TwoPhaseFlowAirTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
   _u_old(valueOld()),
   _T(coupledValue("Temp")),
   _T_old(coupledValue("Temp")),
   _T_id(coupled("Temp")),
   _Pcap(coupledValue("Pcap")),
   _Pcap_old(coupledValue("Pcap")),
   _Pcap_id(coupled("Pcap")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _porosity(getMaterialProperty<Real>("porosity"))
{
}

Real
TwoPhaseFlowAirTimeDerivative::computeQpResidual()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;
  
  // Vapor property
  Real vapor_pressure = 0.001 * R * _T[_qp] / Mw * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure_old = 0.001 * R * _T_old[_qp] / Mw * exp(19.84 - 4976 / _T_old[_qp] + _Pcap_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  // Dry air density value
  Real air_density = Ma / (R * _T[_qp]) * (_u[_qp] - vapor_pressure);
  Real air_density_old = Ma / (R * _T_old[_qp]) * (_u_old[_qp] - vapor_pressure_old);
  
  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _Pcap_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Gas saturation value
  Real Sat_g = 1 - Sat;
  Real Sat_g_old = 1 - Sat_old;

  // Time derivative term
  Real dSat_g_dt = (Sat_g - Sat_g_old) / _dt;
  Real dair_density_dt = (air_density - air_density_old) / _dt;

  return _test[_i][_qp] * _porosity[_qp] * (air_density * dSat_g_dt + Sat_g * dair_density_dt);
}

Real
TwoPhaseFlowAirTimeDerivative::computeQpJacobian()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // Vapor density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] + _Pcap_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  // Vapor pressure value
  Real vapor_pressure = 0.001 * R * _T[_qp] / Mw * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure_old = 0.001 * R * _T_old[_qp] / Mw * exp(19.84 - 4976 / _T_old[_qp] + _Pcap_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  // Dry air density value
  Real air_density = Ma / (R * _T[_qp]) * (_u[_qp] - vapor_pressure);
  Real air_density_old = Ma / (R * _T_old[_qp]) * (_u_old[_qp] - vapor_pressure_old);
  
  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _Pcap_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Gas saturation value
  Real Sat_g = 1 - Sat;
  Real Sat_g_old = 1 - Sat_old;

  // Derivative term
  Real dair_density_dPgas = Ma / (R * _T[_qp]) * _phi[_j][_qp];

  // Time derivative term
  Real dSat_g_dt = (Sat_g - Sat_g_old) / _dt;
  Real dair_density_dt = (air_density - air_density_old) / _dt;

  // Jacobian term
  Real dm = _test[_i][_qp] * _porosity[_qp] * dair_density_dPgas * dSat_g_dt;
      dm += _test[_i][_qp] * _porosity[_qp] * Sat_g * dair_density_dPgas / _dt;

  return dm;
}

Real
TwoPhaseFlowAirTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // Vapor density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_density_old = 0.001 * exp(19.84 - 4976 / _T_old[_qp] + _Pcap_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  Real vapor_pressure = 0.001 * R * _T[_qp] / Mw * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure_old = 0.001 * R * _T_old[_qp] / Mw * exp(19.84 - 4976 / _T_old[_qp] + _Pcap_old[_qp] * Mw / (rho_w * R * _T_old[_qp]));

  // Dry air density value
  Real air_density = Ma / (R * _T[_qp]) * (_u[_qp] - vapor_pressure);
  Real air_density_old = Ma / (R * _T_old[_qp]) * (_u_old[_qp] - vapor_pressure_old);
  
  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real Sat_old = pow((1 + pow((-_alpha[_qp] * _Pcap_old[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Gas saturation value
  Real Sat_g = 1 - Sat;
  Real Sat_g_old = 1 - Sat_old;

  // Time derivative term
  Real dSat_g_dt = (Sat_g - Sat_g_old) / _dt;
  Real dair_density_dt = (air_density - air_density_old) / _dt;

  // Derivative term
  Real dair_density_dPgas = Ma / (R * _T[_qp]) * _phi[_j][_qp];
  Real dvapor_density_dT = (4976 / pow(_T[_qp],2) - _Pcap[_qp] * Mw / (rho_w * R * pow(_T[_qp],2))) * _phi[_j][_qp] * vapor_density;
  Real dvapor_pressure_dT = dvapor_density_dT * R * _T[_qp] / Mw + vapor_density * R / Mw * _phi[_j][_qp];

  Real dair_density_dT = -Ma / (R * _T[_qp]) * dvapor_pressure_dT - (_u[_qp] - vapor_pressure) * Ma / (R * pow(_T[_qp],2)) * _phi[_j][_qp];
  Real dvapor_density_dPcap = Mw / (rho_w * R * _T[_qp]) * vapor_density * _phi[_j][_qp];
  Real dvapor_pressure_dPcap = Mw / (rho_w * R * _T[_qp]) * vapor_pressure * _phi[_j][_qp];

  Real dSat_dPcap = _alpha[_qp] * _lambda[_qp] / (1 - _lambda[_qp]) * pow((-_alpha[_qp] * _Pcap[_qp]),(_lambda[_qp] / (1 - _lambda[_qp]))) * pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp] - 1)) * _phi[_j][_qp];
  Real dSat_g_dPcap = -dSat_dPcap;

  Real dair_density_dPcap = -Ma / (R * _T[_qp]) * dvapor_pressure_dPcap;

  // Jacobian term
  Real dm_T = _test[_i][_qp] * _porosity[_qp] * dair_density_dT * dSat_g_dt;
      dm_T += _test[_i][_qp] * _porosity[_qp] * Sat_g * dair_density_dt;

  Real dm_Pcap = _test[_i][_qp] * _porosity[_qp] * dair_density_dPcap * dSat_g_dt;
      dm_Pcap += _test[_i][_qp] * _porosity[_qp] * air_density * dSat_g_dPcap / _dt;
      dm_Pcap += _test[_i][_qp] * _porosity[_qp] * dSat_g_dPcap * dair_density_dt;
      dm_Pcap += _test[_i][_qp] * _porosity[_qp] * Sat_g * dair_density_dPcap / _dt;

  if (jvar == _T_id)
	return dm_T;
  else if (jvar == _Pcap_id)
	return dm_Pcap;
  else
	return 0.0;  
}

