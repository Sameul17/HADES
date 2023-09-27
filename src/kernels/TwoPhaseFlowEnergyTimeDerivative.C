//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowEnergyTimeDerivative.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include <limits>

registerMooseObject("HadesApp", TwoPhaseFlowEnergyTimeDerivative);

InputParameters
TwoPhaseFlowEnergyTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  params.addCoupledVar("Pcap", -1E7, "The capillary pressure of each quadrature point");
  params.addCoupledVar("Pgas", 1E5, "The gas pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Medium_density", "The density of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Medium_specific_heat", "The specific heat of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The density of water of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_specific_heat", "The specific heat of water of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Gas_specific_heat", "The specific heat of gas of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Porosity", "The porosity of medium of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The lambda of Van Genuchten equation of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The alpha of Van Genuchten equation of each quadrature point");
  

  params.addClassDescription("Derivative of temperature with respect to time. The variable u is temperature in here");
  return params;
}

TwoPhaseFlowEnergyTimeDerivative::TwoPhaseFlowEnergyTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
   _u_old(valueOld()),
   _Pcap(coupledValue("Pcap")),
   _Pcap_old(coupledValue("Pcap")),
   _Pcap_id(coupled("Pcap")),
   _Pgas(coupledValue("Pgas")),
   _Pgas_old(coupledValue("Pgas")),
   _Pgas_id(coupled("Pgas")),
   _Mdensity(getMaterialProperty<Real>("Medium_density")),
   _MCp(getMaterialProperty<Real>("Medium_specific_heat")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _WCp(getMaterialProperty<Real>("Water_specific_heat")),
   _GCp(getMaterialProperty<Real>("Gas_specific_heat")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _porosity(getMaterialProperty<Real>("Porosity"))
{
}

// Variable is Temperature

Real
TwoPhaseFlowEnergyTimeDerivative::computeQpResidual()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  Real vapor_density = 0.001 * exp(19.84 - 4976 / _u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp]));
  Real vapor_pressure = 0.001 * exp(19.84 - 4976 / _u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp])) * R * _u[_qp] / Mw;
  Real air_density = Ma / (R * _u[_qp]) * (_Pgas[_qp] - vapor_pressure);
  Real gas_density = vapor_density + air_density;

  Real rho_Cp_s = _Mdensity[_qp] * _MCp[_qp];
  Real rho_Cp_l = _Wdensity[_qp] * _WCp[_qp];
  Real rho_Cp_g = gas_density * _GCp[_qp];

  Real rho_Cp_eff = (1 - _porosity[_qp]) * rho_Cp_s + _porosity[_qp] * (Sat * rho_Cp_l + (1 - Sat) * rho_Cp_g);
 
  return _test[_i][_qp] * rho_Cp_eff * (_u[_qp] - _u_old[_qp]) / _dt;
}

Real
TwoPhaseFlowEnergyTimeDerivative::computeQpJacobian()
{
  Real R = 8.314; // Universal gas25yy constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  Real vapor_density = 0.001 * exp(19.84 - 4976 / _u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp]));
  Real vapor_pressure = 0.001 * exp(19.84 - 4976 / _u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp])) * R * _u[_qp] / Mw;
  Real air_density = Ma / (R * _u[_qp]) * (_Pgas[_qp] - vapor_pressure);
  Real gas_density = vapor_density + air_density;

  Real rho_Cp_s = _Mdensity[_qp] * _MCp[_qp];
  Real rho_Cp_l = _Wdensity[_qp] * _WCp[_qp];
  Real rho_Cp_g = gas_density * _GCp[_qp];

  Real rho_Cp_eff = (1 - _porosity[_qp]) * rho_Cp_s + _porosity[_qp] * (Sat * rho_Cp_l + (1 - Sat) * rho_Cp_g);

  // Jacobian Derivative term
  Real dvapor_density_dT = (4976 / pow(_u[_qp],2) - _Pcap[_qp] * Mw / (rho_w * R * pow(_u[_qp],2))) * _phi[_j][_qp] * vapor_density;

  Real dvapor_pressure_dT = 0.001 * exp(19.84 - 4976 / _u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp])) * R / Mw * _phi[_j][_qp];
      dvapor_pressure_dT += (4976 / pow(_u[_qp],2) - _Pcap[_qp] * Mw / (rho_w * R * pow(_u[_qp],2))) * _phi[_j][_qp] * vapor_pressure;
  Real dair_density_dT = -Ma / (R * _u[_qp]) * dvapor_pressure_dT - (_Pgas[_qp] - vapor_pressure) * Ma / (R * pow(_u[_qp],2)) * _phi[_j][_qp];
  Real dgas_density_dT = dvapor_density_dT + dair_density_dT;
  Real drho_Cp_eff_dT = _porosity[_qp] * (1 - Sat) * _GCp[_qp] * dgas_density_dT;
 
  Real dm = _test[_i][_qp] * drho_Cp_eff_dT * (_u[_qp] - _u_old[_qp]) / _dt;
      dm += _test[_i][_qp] * rho_Cp_eff * _phi[_j][_qp] / _dt;

  return dm;
}

Real
TwoPhaseFlowEnergyTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Universal gas25yy constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  Real vapor_density = 0.001 * exp(19.84 - 4976/_u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp]));
  Real vapor_pressure = vapor_density * R * _u[_qp] / Mw;
  Real air_density = Ma / (R * _u[_qp]) * (_Pgas[_qp] - vapor_pressure);
  Real gas_density = vapor_density + air_density;

  Real rho_Cp_s = _Mdensity[_qp] * _MCp[_qp];
  Real rho_Cp_l = _Wdensity[_qp] * _WCp[_qp];
  Real rho_Cp_g = gas_density * _GCp[_qp];

  Real rho_Cp_eff = (1 - _porosity[_qp]) * rho_Cp_s + _porosity[_qp] * (Sat * rho_Cp_l + (1 - Sat) * rho_Cp_g);

  // Jacobian derivative term
  Real dvapor_density_dT = (4976 / pow(_u[_qp],2) - _Pcap[_qp] * Mw / (rho_w * R * pow(_u[_qp],2))) * _phi[_j][_qp] * vapor_density;

  Real dvapor_pressure_dT = 0.001 * exp(19.84 - 4976 / _u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp])) * R / Mw * _phi[_j][_qp]; 
      dvapor_pressure_dT += (4976 / pow(_u[_qp],2) - _Pcap[_qp] * Mw / (rho_w * R * pow(_u[_qp],2))) * _phi[_j][_qp] * vapor_pressure;
  Real dair_density_dT = -Ma / (R * _u[_qp]) * dvapor_pressure_dT - (_Pgas[_qp] - vapor_pressure) * Ma / (R * pow(_u[_qp],2)) * _phi[_j][_qp];
  Real dgas_density_dT = dvapor_density_dT + dair_density_dT;
  Real drho_Cp_eff_dT = _porosity[_qp] * (1 - Sat) * _GCp[_qp] * dgas_density_dT;
 
  // OffDiagonalJacobian derivative term
  Real dSat_dPcap = _alpha[_qp] * _lambda[_qp] / (1 - _lambda[_qp]) * pow((-_alpha[_qp] * _Pcap[_qp]),(_lambda[_qp] / (1 - _lambda[_qp]))) * pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp] - 1)) * _phi[_j][_qp];
  Real dSat_g_dPcap = -_alpha[_qp] * _lambda[_qp] / (1 - _lambda[_qp]) * pow((-_alpha[_qp] * _Pcap[_qp]),(_lambda[_qp] / (1 - _lambda[_qp]))) * pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp] - 1)) * _phi[_j][_qp];

  Real dvapor_density_dPcap = Mw / (rho_w * R * _u[_qp]) * _phi[_j][_qp] * vapor_density;
  Real dvapor_pressure_dPcap = Mw / (rho_w * R * _u[_qp]) * _phi[_j][_qp] * vapor_density * R * _u[_qp] / Mw;
  Real dair_density_dPcap = -Ma / (R * _u[_qp]) * dvapor_pressure_dPcap;

  Real dgas_density_dPcap = dvapor_density_dPcap + dair_density_dPcap;
  Real dgas_density_dPgas = Ma / (R * _u[_qp]) * _phi[_j][_qp];

  Real drho_Cp_eff_dPcap = _porosity[_qp] * rho_w * _WCp[_qp] * dSat_dPcap;
      drho_Cp_eff_dPcap += _porosity[_qp] * gas_density * _GCp[_qp] * dSat_g_dPcap;
      drho_Cp_eff_dPcap += _porosity[_qp] * (1 - Sat) * _GCp[_qp] * dgas_density_dPcap;
  Real drho_Cp_eff_dPgas = _porosity[_qp] * (1 - Sat) * _GCp[_qp] * dgas_density_dPgas;

  Real dm_Pcap = _test[_i][_qp] * drho_Cp_eff_dPcap * (_u[_qp] - _u_old[_qp]) / _dt;
  Real dm_Pgas = _test[_i][_qp] * drho_Cp_eff_dPgas * (_u[_qp] - _u_old[_qp]) / _dt;

  if (jvar == _Pcap_id) 
        return dm_Pcap;
  else if (jvar == _Pgas_id)
	return dm_Pgas;
  else
	return 0.0;
}

