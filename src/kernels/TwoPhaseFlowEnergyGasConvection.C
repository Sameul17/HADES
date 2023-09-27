//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowEnergyGasConvection.h"

registerMooseObject("HadesApp", TwoPhaseFlowEnergyGasConvection);

InputParameters
TwoPhaseFlowEnergyGasConvection::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Pcap", -1E4,  "The capillary pressure of each quadrature point");
  params.addCoupledVar("Pgas", 1E5, "The gas pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Intrinsic_permeability", "The intrinsic permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Relative_permeability_gas", "The Relative permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_specific_heat", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Gas_specific_heat", "The specific heat of gas of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Gas_viscosity", "The water viscosity of each quadrature point");
  params.addRequiredParam<RealVectorValue>("gravity", "Gravitational acceleration vector downwards (m/s^2)");

  params.addClassDescription("Calculate energy convection in porous medium");
  return params;
}

TwoPhaseFlowEnergyGasConvection::TwoPhaseFlowEnergyGasConvection(const InputParameters & parameters)
  : Kernel(parameters),
   _Pcap(coupledValue("Pcap")),
   _grad_Pcap(coupledGradient("Pcap")),
   _Pcap_id(coupled("Pcap")),
   _Pgas(coupledValue("Pgas")),
   _grad_Pgas(coupledGradient("Pgas")),
   _Pgas_id(coupled("Pgas")),
   _IPerm(getMaterialProperty<Real>("Intrinsic_permeability")),
   _RPerm(getMaterialProperty<Real>("Relative_permeability_gas")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _WCp(getMaterialProperty<Real>("Water_specific_heat")),
   _GCp(getMaterialProperty<Real>("Gas_specific_heat")),
   _Gviscosity(getMaterialProperty<Real>("Gas_viscosity")),
   _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
TwoPhaseFlowEnergyGasConvection::computeQpResidual()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten

  Real vapor_density = 0.001 * exp(19.84 - 4976/_u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp]));
  Real vapor_pressure = vapor_density * R * _u[_qp] / Mw;
  Real air_density = Ma / (R * _u[_qp]) * (_Pgas[_qp] - vapor_pressure);
  Real gas_density = vapor_density + air_density;

  Real rho_Cp_g = gas_density * _GCp[_qp];

  return _test[_i][_qp] * rho_Cp_g * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_Pgas[_qp] + gas_density * _gravity) * _grad_u[_qp];
}

Real
TwoPhaseFlowEnergyGasConvection::computeQpJacobian()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten

  Real vapor_density = 0.001 * exp(19.84 - 4976/_u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp]));
  Real vapor_pressure = vapor_density * R * _u[_qp] / Mw;
  Real air_density = Ma / (R * _u[_qp]) * (_Pgas[_qp] - vapor_pressure);
  Real gas_density = vapor_density + air_density;

  Real rho_Cp_g = gas_density * _GCp[_qp];

  return _test[_i][_qp] * rho_Cp_g * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_Pgas[_qp] + gas_density * _gravity) * _grad_phi[_j][_qp];
}

Real
TwoPhaseFlowEnergyGasConvection::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten

  Real vapor_density = 0.001 * exp(19.84 - 4976/_u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp]));
  Real vapor_pressure = vapor_density * R * _u[_qp] / Mw;
  Real air_density = Ma / (R * _u[_qp]) * (_Pgas[_qp] - vapor_pressure);
  Real gas_density = vapor_density + air_density;

  Real rho_Cp_g = gas_density * _GCp[_qp];

  // Derivative term
  Real dvapor_density_dPcap = 0.001 * Mw / (_Wdensity[_qp] * R * _u[_qp]) * exp(19.84 - 4976 / _u[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _u[_qp])) * _phi[_j][_qp];
  Real dair_density_dPcap = -Ma / Mw * dvapor_density_dPcap;

  Real dgas_density_dPcap = dvapor_density_dPcap + dair_density_dPcap;
  Real dgas_density_dPgas = Ma / (R * _u[_qp]) * _phi[_j][_qp];

  Real drho_Cp_g_dPcap = dgas_density_dPcap * _GCp[_qp];
  Real drho_Cp_g_dPgas = dgas_density_dPgas * _GCp[_qp];
  
  Real dm_Pcap = _test[_i][_qp] * drho_Cp_g_dPcap * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_Pgas[_qp] + gas_density * _gravity) * _grad_u[_qp];
      dm_Pcap += _test[_i][_qp] * rho_Cp_g * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (dgas_density_dPcap * _gravity) * _grad_u[_qp];

  Real dm_Pgas = _test[_i][_qp] * drho_Cp_g_dPgas * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_Pgas[_qp] + gas_density * _gravity) * _grad_u[_qp];
      dm_Pgas += _test[_i][_qp] * rho_Cp_g * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_phi[_j][_qp] + dgas_density_dPgas * _gravity) * _grad_u[_qp];

  if (jvar == _Pcap_id)
    return dm_Pcap;
  else if (jvar == _Pgas_id)
    return dm_Pgas;
  else
	  return 0.0;
}
