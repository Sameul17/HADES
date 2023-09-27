//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowAirAdvectiveFluxTHMC.h"

registerMooseObject("HadesApp", TwoPhaseFlowAirAdvectiveFluxTHMC);

InputParameters
TwoPhaseFlowAirAdvectiveFluxTHMC::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point");
  params.addCoupledVar("Pcap", -1E4, "The gas pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Gas_viscosity", "The dynamic viscosity of water of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The Van Genuchten parameter of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The Van Genuchten parameter of each quadrature point");
  params.addRequiredParam<RealVectorValue>("gravity","Gravitational acceleration vector downwards (m/s^2)");

  params.addClassDescription("Calculate water advection in porous medium");
  return params;
}

TwoPhaseFlowAirAdvectiveFluxTHMC::TwoPhaseFlowAirAdvectiveFluxTHMC(const InputParameters & parameters)
  : Kernel(parameters),
   _T(coupledValue("Temp")),
   _grad_T(coupledGradient("Temp")),
   _T_id(coupled("Temp")),
   _Pcap(coupledValue("Pcap")),
   _grad_Pcap(coupledGradient("Pcap")),
   _Pcap_id(coupled("Pcap")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _Gviscosity(getMaterialProperty<Real>("Gas_viscosity")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
   _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
TwoPhaseFlowAirAdvectiveFluxTHMC::computeQpResidual()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Vapor density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure = vapor_density * R * _T[_qp] / Mw;
  
  // Dry air density value
  Real air_density = Ma / (R * _T[_qp]) * (_u[_qp] - vapor_pressure);
  Real air_pressure = air_density * R * _T[_qp] / Ma;
  Real gas_density = vapor_density + air_density;

  // Permeability value
  Real kk_rel = 5.695E-11 * pow(1.44 * (1 - Sat) ,4.3);

//  return _grad_test[_i][_qp] * air_density * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
  return _grad_test[_i][_qp] * air_density * kk_rel / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
}

Real
TwoPhaseFlowAirAdvectiveFluxTHMC::computeQpJacobian()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Vapor density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure = vapor_density * R * _T[_qp] / Mw;
  
  // Dry air density value
  Real air_density = Ma / (R * _T[_qp]) * (_u[_qp] - vapor_pressure);
  Real air_pressure = air_density * R * _T[_qp] / Ma;
  Real gas_density = vapor_density + air_density;

  // Derivative term
  Real dair_density_dPgas = Ma / (R * _T[_qp]) * _phi[_j][_qp];
  Real dgas_density_dPgas = Ma / (R * _T[_qp]) * _phi[_j][_qp];

  // Permeability value
  Real kk_rel = 5.695E-11 * pow(1.44 * (1 - Sat) ,4.3);

  // Jacobian term
//  Real dm = _grad_test[_i][_qp] * dair_density_dPgas * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
//  dm += _grad_test[_i][_qp] * air_density * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_phi[_j][_qp] - dgas_density_dPgas * _gravity);
  Real dm = _grad_test[_i][_qp] * dair_density_dPgas * kk_rel / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
  dm += _grad_test[_i][_qp] * air_density * kk_rel / _Gviscosity[_qp] * (_grad_phi[_j][_qp] - dgas_density_dPgas * _gravity);

  return dm;
}

Real
TwoPhaseFlowAirAdvectiveFluxTHMC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Vapor density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp]) * exp(_Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure = vapor_density * R * _T[_qp] / Mw;
  
  // Dry air density value
  Real air_density = Ma / (R * _T[_qp]) * (_u[_qp] - vapor_pressure);
  Real air_pressure = air_density * R * _T[_qp] / Ma;
  Real gas_density = vapor_density + air_density;

  // Permeability value
  Real kk_rel = 5.695E-11 * pow(1.44 * (1 - Sat) ,4.3);

  // Derivative term
  Real dair_density_dPgas = Ma / (R * _T[_qp]) * _phi[_j][_qp];
  Real dgas_density_dPgas = Ma / (R * _T[_qp]) * _phi[_j][_qp];

  Real dvapor_density_dT = ((4976 / (_T[_qp] * _T[_qp])) + Mw / (_Wdensity[_qp] * R * _T[_qp] *_T[_qp])) * vapor_density * _phi[_j][_qp];
  Real dvapor_pressure_dT = vapor_density * R / Mw * _phi[_j][_qp] + R * _T[_qp] / Mw * dvapor_density_dT;

  Real dvapor_density_dPcap = Mw / (_Wdensity[_qp] * R * _T[_qp]) * vapor_density;
  Real dvapor_pressure_dPcap = dvapor_density_dPcap * R * _T[_qp] / Mw;

  Real dair_density_dT = -Ma / (R * pow(_T[_qp],2)) * (_u[_qp] - vapor_pressure) * _phi[_j][_qp] - Ma / (R * _T[_qp]) * dvapor_density_dT;
  Real dgas_density_dT = dair_density_dT + dvapor_density_dT;

  Real dair_density_dPcap = -Ma / (R * _T[_qp]) * dvapor_pressure_dPcap;
  Real dgas_density_dPcap = dair_density_dPcap + dvapor_density_dPcap;

  // Jacobian term
//  Real dm_T = _grad_test[_i][_qp] * dair_density_dT * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
//      dm_T += _grad_test[_i][_qp] * air_density * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * -dgas_density_dT * _gravity;

//  Real dm_Pcap = _grad_test[_i][_qp] * dair_density_dPcap * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
//      dm_Pcap += _grad_test[_i][_qp] * air_density * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * -dgas_density_dPcap * _gravity;

  Real dm_T = _grad_test[_i][_qp] * dair_density_dT * kk_rel / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
      dm_T += _grad_test[_i][_qp] * air_density * kk_rel / _Gviscosity[_qp] * -dgas_density_dT * _gravity;

  Real dm_Pcap = _grad_test[_i][_qp] * dair_density_dPcap * kk_rel / _Gviscosity[_qp] * (_grad_u[_qp] - gas_density * _gravity);
      dm_Pcap += _grad_test[_i][_qp] * air_density * kk_rel / _Gviscosity[_qp] * -dgas_density_dPcap * _gravity;

  if (jvar == _T_id)
       return dm_T;
  else if (jvar == _Pcap_id)
       return dm_Pcap; 
  else
       return 0.0;
}
