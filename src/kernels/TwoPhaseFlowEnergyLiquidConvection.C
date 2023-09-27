//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowEnergyLiquidConvection.h"

registerMooseObject("HadesApp", TwoPhaseFlowEnergyLiquidConvection);

InputParameters
TwoPhaseFlowEnergyLiquidConvection::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Pcap", -1E4, "The capillary pressure of each quadrature point");
  params.addCoupledVar("Pgas", 1E5, "The gas pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Intrinsic_permeability", "The intrinsic permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Relative_permeability", "The Relative permeability of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_specific_heat", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_viscosity", "The water viscosity of each quadrature point");
  params.addRequiredParam<RealVectorValue>("gravity", "Gravitational acceleration vector downwards (m/s^2)");

  params.addClassDescription("Calculate energy convection in porous medium");
  return params;
}

TwoPhaseFlowEnergyLiquidConvection::TwoPhaseFlowEnergyLiquidConvection(const InputParameters & parameters)
  : Kernel(parameters),
   _Pcap(coupledValue("Pcap")),
   _grad_Pcap(coupledGradient("Pcap")),
   _Pcap_id(coupled("Pcap")),
   _Pgas(coupledValue("Pgas")),
   _grad_Pgas(coupledGradient("Pgas")),
   _Pgas_id(coupled("Pgas")),
   _IPerm(getMaterialProperty<Real>("Intrinsic_permeability")),
   _RPerm(getMaterialProperty<Real>("Relative_permeability")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _WCp(getMaterialProperty<Real>("Water_specific_heat")),
   _Wviscosity(getMaterialProperty<Real>("Water_viscosity")),
   _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
TwoPhaseFlowEnergyLiquidConvection::computeQpResidual()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;
  // The equaiton of relation between degree of saturaiton and capillary pressure given by Van Genuchten

  Real vapor_density = 0.001 * exp(19.84 - 4976/_u[_qp] + _Pcap[_qp] * Mw / (_Wdensity[_qp] * R * _u[_qp]));
  Real vapor_pressure = vapor_density * R * _u[_qp] / Mw;
  Real air_density = Ma / (R * _u[_qp]) * (_Pgas[_qp] - vapor_pressure);
  Real gas_density = vapor_density + air_density;

  Real rho_Cp_l = rho_w * _WCp[_qp];

  return _test[_i][_qp] * rho_Cp_l * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_Pgas[_qp] - _grad_Pcap[_qp] + rho_w * _gravity) * _grad_u[_qp];
}

Real
TwoPhaseFlowEnergyLiquidConvection::computeQpJacobian()
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

  Real rho_Cp_l = rho_w * _WCp[_qp];

  return _test[_i][_qp] * rho_Cp_l * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * (_grad_Pgas[_qp] - _grad_Pcap[_qp] + rho_w * _gravity) * _grad_phi[_j][_qp];
}

Real
TwoPhaseFlowEnergyLiquidConvection::computeQpOffDiagJacobian(unsigned int jvar)
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

  Real rho_Cp_l = rho_w * _WCp[_qp];
  
  Real dm_Pcap = -_test[_i][_qp] * rho_Cp_l * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * _grad_phi[_j][_qp] * _grad_u[_qp];
  Real dm_Pgas = _test[_i][_qp] * rho_Cp_l * _IPerm[_qp] * _RPerm[_qp] / _Wviscosity[_qp] * _grad_phi[_j][_qp] * _grad_u[_qp];

  if (jvar == _Pcap_id)
    return dm_Pcap;
  else if (jvar == _Pgas_id)
    return dm_Pgas;
  else
	  return 0.0;
}
