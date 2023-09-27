//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowVaporAdvectiveFluxTHMC.h"

registerMooseObject("HadesApp", TwoPhaseFlowVaporAdvectiveFluxTHMC);

InputParameters
TwoPhaseFlowVaporAdvectiveFluxTHMC::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point");
  params.addCoupledVar("Pgas", 1E5, "The water pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The ags density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Gas_viscosity", "The dynamic viscosity of gas of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The Van Genuchten parameter of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The Van Genuchten parameter of each quadrature point");
  params.addRequiredParam<RealVectorValue>("gravity","Gravitational acceleration vector downwards (m/s^2)");

  params.addClassDescription("Calculate gas advection in porous medium");
  return params;
}

TwoPhaseFlowVaporAdvectiveFluxTHMC::TwoPhaseFlowVaporAdvectiveFluxTHMC(const InputParameters & parameters)
  : Kernel(parameters),
   _T(coupledValue("Temp")),
   _T_id(coupled("Temp")),
   _Pgas(coupledValue("Pgas")),
   _grad_Pgas(coupledGradient("Pgas")),
   _Pgas_id(coupled("Pgas")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _Gviscosity(getMaterialProperty<Real>("Gas_viscosity")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
   _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
TwoPhaseFlowVaporAdvectiveFluxTHMC::computeQpResidual()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp])) * R * _T[_qp] / Mw;
  Real air_density = (_Pgas[_qp] - vapor_pressure) * Ma / (R * _T[_qp]);
  Real gas_density = vapor_density + air_density;

  Real kk_rel = 5.695E-11 * pow(1.44 * (1 - Sat),4.3);

//  return _grad_test[_i][_qp] * vapor_density * _IPerm[_qp] * _RPerm[_qp] / _Gviscosity[_qp] * (_grad_Pgas[_qp] - gas_density * _gravity);
  return _grad_test[_i][_qp] * vapor_density * kk_rel / _Gviscosity[_qp] * (_grad_Pgas[_qp] - gas_density * _gravity);
}

Real
TwoPhaseFlowVaporAdvectiveFluxTHMC::computeQpJacobian()
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;

  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp])) * R * _T[_qp] / Mw;
  Real air_density = (_Pgas[_qp] - vapor_pressure) * Ma / (R * _T[_qp]);
  Real gas_density = vapor_density + air_density;

  // Primary Variable Derivative term
  Real dvapor_density_dPcap = Mw / (_Wdensity[_qp] * R * _T[_qp]) * vapor_density * _phi[_j][_qp];
  Real dvapor_pressure_dPcap = Mw / (_Wdensity[_qp] * R * _T[_qp]) * vapor_pressure * _phi[_j][_qp];
  Real dair_density_dPcap = -dvapor_pressure_dPcap * Ma / (R * _T[_qp]);
  Real dgas_density_dPcap = dvapor_density_dPcap + dair_density_dPcap; 
  Real kk_rel = 5.695E-11 * pow(1.44 * (1 - Sat),4.3);

// Jacobian term
  Real dm = _grad_test[_i][_qp] * dvapor_density_dPcap * kk_rel / _Gviscosity[_qp] * (_grad_Pgas[_qp] - gas_density * _gravity);
      dm -= _grad_test[_i][_qp] * vapor_density * kk_rel / _Gviscosity[_qp] * dgas_density_dPcap * _gravity;

  return dm; 
}

Real
TwoPhaseFlowVaporAdvectiveFluxTHMC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Universal gas constant [J/(kgK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Ma = 0.0289647; // Molar mass of dry air [kg/mol]
  Real rho_w = 1E3;
  Real T2 = pow(_T[_qp],2); // Temperature square [K^2]

  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real vapor_pressure = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp])) * R * _T[_qp] / Mw;
  Real air_density = (_Pgas[_qp] - vapor_pressure) * Ma / (R * _T[_qp]);
  Real gas_density = vapor_density + air_density;

  // Primary Variable Derivative term
  Real dvapor_density_dPcap = Mw / (rho_w * R * _T[_qp]) * vapor_density * _phi[_j][_qp];
  Real dvapor_pressure_dPcap = Mw / (rho_w * R * _T[_qp]) * vapor_pressure * _phi[_j][_qp];
  Real dair_density_dPcap = -dvapor_pressure_dPcap * Ma / (R * _T[_qp]);
  Real dgas_density_dPcap = dvapor_density_dPcap + dair_density_dPcap; 

  // Secondary Variable Derivative term
  Real dvapor_density_dT = (4976 / T2 + _u[_qp] * Mw / (rho_w * R * T2)) * vapor_density * _phi[_j][_qp];
  Real dvapor_pressure_dT = (4976 / T2 + _u[_qp] * Mw / (rho_w * R * T2)) * vapor_pressure * R * _T[_qp] / Mw * _phi[_j][_qp] ;
      dvapor_pressure_dT += vapor_pressure * R / Mw * _phi[_j][_qp];
  Real dair_density_dT = -dvapor_pressure_dT * Ma / (R * _T[_qp]) + (_Pgas[_qp] - vapor_pressure) * Ma / (R * pow(_T[_qp],2)) * _phi[_j][_qp];
  Real dgas_density_dT = dvapor_density_dT + dair_density_dT;

  Real dvapor_density_dPgas = 0;
  Real dair_density_dPgas = Ma / (R * _T[_qp]) * _phi[_j][_qp];
  Real dgas_density_dPgas = dvapor_density_dPgas + dair_density_dPgas;
  Real kk_rel = 5.695E-11 * pow(1.44 * (1 - Sat),4.3);


  Real dm_T = _grad_test[_i][_qp] * dvapor_density_dT * kk_rel / _Gviscosity[_qp] * (_grad_Pgas[_qp] - gas_density * _gravity);
      dm_T -= _grad_test[_i][_qp] * vapor_density * kk_rel / _Gviscosity[_qp] * (dgas_density_dT * _gravity);

  Real dm_Pgas = _grad_test[_i][_qp] * vapor_density * kk_rel / _Gviscosity[_qp] * (_grad_phi[_j][_qp] - dgas_density_dPgas * _gravity);


  if (jvar == _T_id)
	return dm_T;
  else if (jvar == _Pgas_id)
	return dm_Pgas;
  else	  	
	return 0.0;
}
