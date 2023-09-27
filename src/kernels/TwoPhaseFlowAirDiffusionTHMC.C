//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseFlowAirDiffusionTHMC.h"

registerMooseObject("HadesApp", TwoPhaseFlowAirDiffusionTHMC);

InputParameters
TwoPhaseFlowAirDiffusionTHMC::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addCoupledVar("Temp", 300, "The temperature of each quadrature point");
  params.addCoupledVar("Pcap", -1E5, "The gas pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Water_density", "The water density of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("porosity", "The Bentonite porosity of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The Van Genuchten parameter of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The Van Genuchten parameter of each quadrature point");
  params.addClassDescription("Calculate the Gas Diffusion behavior in porous medium");
  return params;
}

TwoPhaseFlowAirDiffusionTHMC::TwoPhaseFlowAirDiffusionTHMC(const InputParameters & parameters)
  : Kernel(parameters),
   _T(coupledValue("Temp")),
   _grad_T(coupledGradient("Temp")),
   _T_id(coupled("Temp")),
   _Pcap(coupledValue("Pcap")),
   _grad_Pcap(coupledGradient("Pcap")),
   _Pcap_id(coupled("Pcap")),
   _Wdensity(getMaterialProperty<Real>("Water_density")),
   _porosity(getMaterialProperty<Real>("porosity")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda"))
{
}

Real
TwoPhaseFlowAirDiffusionTHMC::computeQpResidual()
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor
  Real rho_w = 1E3;

  // Saturation value
  Real Sat = pow(1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp]))),(-_lambda[_qp]));

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] / (rho_w * Rv * _T[_qp]));
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // The saturated water vapor density

  Real Theta = exp((_Pcap[_qp] + _u[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
//  Real Theta = exp((_Pcap[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
//  Real drho_dT = 4976 / _T[_qp] * rho_vS; // The derivative term of drho_vs / dtemp
  Real drho_dT = 4976 / pow(_T[_qp],2) * rho_vS; // The derivative term of drho_vs / dtemp
  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)

  // Vapor Diffusion related coefficient
  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600 / 10000; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 5.9E-6 * pow(_T[_qp] - 273.15, 2.3) / _u[_qp] * 365*24*3600;


  Real D_pv = D_v * vapor_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity
  Real D_Tv = D_v * (Theta * drho_dT - vapor_density * (_Pcap[_qp] + _u[_qp]) / (rho_w * Rv * _T[_qp] * _T[_qp])); // Thermal vapor diffusivity
  
  return -_grad_test[_i][_qp] * (D_pv * (_grad_Pcap[_qp] - _grad_u[_qp]) + f_Tv * D_Tv * _grad_T[_qp]);
}

Real
TwoPhaseFlowAirDiffusionTHMC::computeQpJacobian()
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor
  Real rho_w = 1E3;

  // Saturation value
  Real Sat = pow(1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp]))),(-_lambda[_qp]));

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // The saturated water vapor density

  Real Theta = exp((_Pcap[_qp] + _u[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
//  Real Theta = exp((_Pcap[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
//  Real drho_dT = 4976 / _T[_qp] * rho_vS; // The derivative term of drho_vs / dtemp
  Real drho_dT = 4976 / pow(_T[_qp],2) * rho_vS; // The derivative term of drho_vs / dtemp
  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)
  
  // Vapor Diffusion related coefficient
  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600 / 10000; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 5.9E-6 * pow(_T[_qp] - 273.15, 2.3) / _u[_qp] * 365*24*3600;

  Real D_pv = D_v * vapor_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity
//  Real D_Tv = D_v * (Theta * drho_dT - vapor_density * (_Pcap[_qp] - _u[_qp]) / (_Wdensity[_qp] * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity

  // Derivative term
//  Real dTheta_dPgas = -1 / (_Wdensity[_qp] * Rv * _T[_qp]) * Theta * _phi[_j][_qp];
//  Real dD_Tv_dPgas = D_v * (dTheta_dPgas * drho_dT + vapor_density / (_Wdensity[_qp] * Rv * pow(_T[_qp],2)) * _phi[_j][_qp]);
  
  //Diagonal Jacobian

  return _grad_test[_i][_qp] * D_pv * _grad_phi[_j][_qp];
}

Real
TwoPhaseFlowAirDiffusionTHMC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314; // Gas cosntant [J/(molK)]
  Real Mw = 0.01801528; // Molar mass of water [kg/mol]
  Real Rv = 461.5; // The specific gas constant [J/(kgK)]
  Real vT = 0.8; //v is the mass flow factor, T is the tortuosity factor
  Real rho_w = 1E3;

  // Saturation value
  Real Sat = pow(1 + pow((-_alpha[_qp] * _Pcap[_qp]),(1 / (1 - _lambda[_qp]))),(-_lambda[_qp]));

  // Gas density value
  Real vapor_density = 0.001 * exp(19.84 - 4976 / _T[_qp] + _Pcap[_qp] * Mw / (rho_w * R * _T[_qp]));
  Real rho_vS = 0.001 * exp(19.84 - 4976 / _T[_qp]); // The saturated water vapor density

  Real Theta = exp((_Pcap[_qp] + _u[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
//  Real Theta = exp((_Pcap[_qp]) / (rho_w * Rv * _T[_qp])); // The relative humidity
//  Real drho_dT = 4976 / _T[_qp] * rho_vS; // The derivative term of drho_vs / dtemp
  Real drho_dT = 4976 / pow(_T[_qp],2) * rho_vS; // The derivative term of drho_vs / dtemp
  Real f_Tv = 1.06; //The thermal diffusion enhancement factor for the thermal diffusion term (f_Tv = 1 + T / 4976)
  
  // Vapor Diffusion related coefficient
//  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600 / 10000; // Vapor Diffusion Coefficient [m2/yr]
  Real D_v = 2.16E-5 * vT * _porosity[_qp] * (1 - Sat) * pow((_T[_qp] / 273),1.8) * 365 * 24 * 3600; // Vapor Diffusion Coefficient [m2/yr]
//  Real D_v = 5.9E-6 * pow(_T[_qp] - 273.15, 2.3) / _u[_qp] * 365*24*3600;

  Real D_pv = D_v * vapor_density / (rho_w * Rv * _T[_qp]); // The isothermal vapor diffusivity
  Real D_Tv = D_v * (Theta * drho_dT - vapor_density * (_Pcap[_qp] + _u[_qp]) / (rho_w * Rv * pow(_T[_qp],2))); // Thermal vapor diffusivity

   Real dm_T = -_grad_test[_i][_qp](0) * f_Tv * D_Tv * _grad_phi[_j][_qp](0);
   Real dm_Pcap = -_grad_test[_i][_qp](0) * D_pv * _grad_phi[_j][_qp](0);

//  Real dm_T = -_grad_test[_i][_qp] * dD_pv_dT * (_grad_u[_qp] - _grad_Pcap[_qp]);
//  dm_T += -_grad_test[_i][_qp] * f_Tv * dD_Tv_dT * _grad_T[_qp];
//  dm_T += -_grad_test[_i][_qp] * f_Tv * D_Tv * _grad_phi[_j][_qp];

//  Real dm_Pcap = -_grad_test[_i][_qp] * dD_pv_dPcap * (_grad_u[_qp] - _grad_Pcap[_qp]);
//  dm_Pcap += _grad_test[_i][_qp] * D_pv * _grad_phi[_j][_qp];
//  dm_Pcap -= _grad_test[_i][_qp] * f_Tv * dD_Tv_dPcap * _grad_T[_qp];

  if (jvar == _T_id)
	  return dm_T;
  else if (jvar == _Pcap_id)
	  return dm_Pcap; 
  else
	  return 0.0;
}
