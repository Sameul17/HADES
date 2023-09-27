//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowVaporDiffusionTHMC.h"

registerMooseObject("HadesApp", TwoPhaseFlowVaporDiffusionTHMC);

InputParameters
TwoPhaseFlowVaporDiffusionTHMC::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point");
  params.addCoupledVar("Pgas", 1E5, "The gas pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("porosity", "The Bentonite porosity of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The Van Genuchten parameter of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The Van Genuchten parameter of each quadrature point");
  params.addClassDescription("Calculate the Gas Diffusion behavior in porous medium");
  return params;
}

TwoPhaseFlowVaporDiffusionTHMC::TwoPhaseFlowVaporDiffusionTHMC(const InputParameters & parameters)
  : Kernel(parameters),
   _T(coupledValue("Temp")),
   _grad_T(coupledGradient("Temp")),
   _T_id(coupled("Temp")),
   _Pgas(coupledValue("Pgas")),
   _grad_Pgas(coupledGradient("Pgas")),
   _Pgas_id(coupled("Pgas")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _porosity(getMaterialProperty<Real>("porosity")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda"))
{
}

Real
TwoPhaseFlowVaporDiffusionTHMC::computeQpResidual()
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor
  Real rho_w = 1E3; // Water density [kg/m3]

  // Saturation value
  Real Sat = pow((1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 5.9E-6 * pow(_T[_qp] - 273.15, 2.3) / _Pgas[_qp] * 365*24*3600;

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real D_pv = D_v * vapor_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity

  Real Theta = exp((_u[_qp] + _Pgas[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // The saturated water vapor density
  Real drho_dT = 4976 / pow(_T[_qp],2) * rho_vS; // The derivative term of drho_vS / dtemp
  Real D_Tv = D_v * (Theta * drho_dT - vapor_density * (_u[_qp] + _Pgas[_qp]) / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity

  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)

  return _grad_test[_i][_qp] * (D_pv * (_grad_u[_qp] - _grad_Pgas[_qp]) + f_Tv * D_Tv * _grad_T[_qp]);
}

Real
TwoPhaseFlowVaporDiffusionTHMC::computeQpJacobian()
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor
  Real rho_w = 1E3; // Water density [kg/m3]

  // Saturation value
  Real Sat = pow(1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp]))),(-_lambda[_qp]));
  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 5.9E-6 * pow(_T[_qp] - 273.15, 2.3) / _Pgas[_qp] * 365*24*3600;

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real D_pv = D_v * vapor_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity

  Real Theta = exp((_u[_qp] + _Pgas[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // The saturated water vapor density
  Real drho_dT = 4976 / pow(_T[_qp],2) * rho_vS; // The derivative term of drho_vS / dtemp
  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)

  // Vapor Diffusion related coefficient
  Real D_Tv = D_v * (Theta * drho_dT - vapor_density * (_u[_qp] + _Pgas[_qp]) / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity
//  Real D_Tv = D_v * (Theta * drho_dT - vapor_density * _u[_qp] / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity

  // Derivative term
//  Real dvapor_density_dPcap = Mw / (_Wdensity[_qp] * R * _T[_qp]) * vapor_density * _phi[_j][_qp];
//  Real dD_pv_dPcap = D_v / (_Wdensity[_qp] * R * _T[_qp]) * dvapor_density_dPcap;

//  Real dTheta_dPcap = 1 / (_Wdensity[_qp] * Rv * _T[_qp]) * Theta * _phi[_j][_qp]; // dTheta / dPcap
//  Real dD_Tv_dPcap = D_v * (drho_dT * dTheta_dPcap - (_Pgas[_qp] - _u[_qp]) / (_Wdensity[_qp] * Rv * pow(_T[_qp],2)) * dvapor_density_dPcap - vapor_density / (_Wdensity[_qp] * Rv * pow(_T[_qp],2)) * -_phi[_j][_qp]); // dD_Tv / dPcap
  
  //Diagonal Jacobian
//  Real dm = _grad_test[_i][_qp] * dD_pv_dPcap * (_grad_Pgas[_qp] - _grad_u[_qp]);
//       dm += _grad_test[_i][_qp] * D_pv * -_grad_phi[_j][_qp];
//       dm += _grad_test[_i][_qp] * dD_Tv_dPcap * f_Tv * _grad_T[_qp];

//  return dm;
  return _grad_test[_i][_qp] * D_pv * _grad_phi[_j][_qp];
}

Real
TwoPhaseFlowVaporDiffusionTHMC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor
  Real rho_w = 1E3; // Water density [kg/m3]

  // Saturation value
  Real Sat = pow(1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp]))),(-_lambda[_qp]));

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // The saturated water vapor density

  Real Theta = exp((_u[_qp] + _Pgas[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
  Real drho_dT = 4976 / pow(_T[_qp],2) * rho_vS; // The derivative term of drho_vS / dtemp
  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)

  // Vapor Diffusion related coefficient
  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 5.9E-6 * pow(_T[_qp] - 273.15, 2.3) / _Pgas[_qp] * 365*24*3600;

  Real D_pv = D_v * vapor_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity
  Real D_Tv = D_v * (Theta * drho_dT - vapor_density * (_u[_qp] + _Pgas[_qp]) / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity

  Real dm_T = _grad_test[_i][_qp](0) * f_Tv * D_Tv * _grad_phi[_j][_qp](0);
  Real dm_Pgas = -_grad_test[_i][_qp](0) * D_pv * _grad_phi[_j][_qp](0);

  if (jvar == _T_id)
	  return dm_T;
  else if (jvar == _Pgas_id)
	  return dm_Pgas; 
  else
	  return 0.0;
}
