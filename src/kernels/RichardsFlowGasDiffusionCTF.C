//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RichardsFlowGasDiffusionCTF.h"

registerMooseObject("HadesApp", RichardsFlowGasDiffusionCTF);

InputParameters
RichardsFlowGasDiffusionCTF::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("porosity", "The Bentonite porosity of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The Van Genuchten parameter of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The Van Genuchten parameter of each quadrature point");
  params.addClassDescription("Calculate the Gas Diffusion behavior in porous medium");
  return params;
}

RichardsFlowGasDiffusionCTF::RichardsFlowGasDiffusionCTF(const InputParameters & parameters)
  : Kernel(parameters),
   _T(coupledValue("Temp")),
   _grad_T(coupledGradient("Temp")),
   _T_id(coupled("Temp")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _porosity(getMaterialProperty<Real>("porosity")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda"))
{
}

Real
RichardsFlowGasDiffusionCTF::computeQpResidual()
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real rho_w = 1E3; // Water density [kg/m3]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor

  Real Sat = pow(1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp]))),-_lambda[_qp]);
  Real D_v = 5.9E-12 * pow((_T[_qp] - 273.15) / 0.1,2.3) * 365 * 24 * 3600; // The effective vapor diffusion coefficient
//  Real D_v = 5.9E-12 * pow(_T[_qp] - 273.15,2.3) / 0.1 * 365 * 24 * 3600; // The effective vapor diffusion coefficient
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600 / 10000; // The effective vapor diffusion coefficient [m2/yr]
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // The effective vapor diffusion coefficient [m2/yr]
  
  Real gas_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real D_pv = D_v * gas_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity

  Real Theta = exp(_u[_qp] / (rho_w * Rv * _T[_qp])); // The relative humidity
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // Saturated vapor density
  Real drho_vS_dT = 4976 / pow(_T[_qp],2) * rho_vS;
//  Real drho_vS_dT = 4976 / _T[_qp] * rho_vS;
  Real D_Tv = D_v * (Theta * drho_vS_dT - gas_density * _u[_qp] / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity

  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)

//  return _grad_test[_i][_qp] * gas_diffusion;
  return _grad_test[_i][_qp] * (D_pv * _grad_u[_qp] + f_Tv * D_Tv * _grad_T[_qp]);
}

Real
RichardsFlowGasDiffusionCTF::computeQpJacobian()
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real rho_w = 1E3; // Water density [g/m3]
//  Real D_v = 763.2; // The effective vapor diffusion coefficient [m2/yr]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor

  Real Sat = pow(1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp]))),(-_lambda[_qp]));
  Real D_v = 5.9E-12 * pow((_T[_qp] - 273.15) / 0.1,2.3) * 365 * 24 * 3600; // The effective vapor diffusion coefficient
//  Real D_v = 5.9E-12 * pow(_T[_qp] - 273.15,2.3) / 0.1 * 365 * 24 * 3600; // The effective vapor diffusion coefficient
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600 / 10000; // The effective vapor diffusion coefficient [m2/yr]
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // The effective vapor diffusion coefficient [m2/yr]
  
  Real gas_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real D_pv = D_v * gas_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity

  Real Theta = exp(_u[_qp] / (rho_w * Rv * _T[_qp])); // The relative humidity
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // Saturated vapor density
  Real drho_vS_dT = 4976 / pow(_T[_qp],2) * rho_vS;
//  Real drho_vS_dT = 4976 / _T[_qp] * rho_vS;
  Real D_Tv = D_v * (Theta * drho_vS_dT - gas_density * _u[_qp] / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity
  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)

//  Real dgas_diffusion_dVar = D_pv * _grad_phi[_j][_qp](0);
      
  return _grad_test[_i][_qp] * D_pv * _grad_phi[_j][_qp];
}

Real
RichardsFlowGasDiffusionCTF::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real rho_w = 1E3; // Water density [g/m3]
//  Real D_v = 763.2; // The effective vapor diffusion coefficient [m2/yr]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor
  Real Sat = pow(1 + pow((-_alpha[_qp] * _u[_qp]),(1 / (1 - _lambda[_qp]))),(-_lambda[_qp]));
  Real D_v = 5.9E-12 * pow((_T[_qp] - 273.15) / 0.1,2.3) * 365 * 24 * 3600; // The effective vapor diffusion coefficient
//  Real D_v = 5.9E-12 * pow(_T[_qp] - 273.15,2.3) / 0.1 * 365 * 24 * 3600; // The effective vapor diffusion coefficient
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600 / 10000; // The effective vapor diffusion coefficient [m2/yr]
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // The effective vapor diffusion coefficient [m2/yr]
  
  Real gas_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _u[_qp] / (rho_w * Rv * _T[_qp]));
  Real D_pv = D_v * gas_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity

  Real Theta = exp(_u[_qp] / (rho_w * Rv * _T[_qp])); // The relative humidity
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // Saturated vapor density
  Real drho_vS_dT = 4976 / pow(_T[_qp],2) * rho_vS;
//  Real drho_vS_dT = 4976 / _T[_qp] * rho_vS;
  Real D_Tv = D_v * (Theta * drho_vS_dT - gas_density * _u[_qp] / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity
  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)

  Real dgas_diffusion_dT = f_Tv * D_Tv * _grad_phi[_j][_qp](0);

  if (jvar == _T_id)
	  return _grad_test[_i][_qp](0) * dgas_diffusion_dT;
  else
	  return 0.0;
}
