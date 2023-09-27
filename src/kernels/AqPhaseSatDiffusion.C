//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AqPhaseSatDiffusion.h"

#include "MooseVariable.h"

#include "libmesh/quadrature.h"

#include <limits>

registerMooseObject("HadesApp", AqPhaseSatDiffusion);

InputParameters
AqPhaseSatDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addCoupledVar("T", 300, "The temperature of each quadrature point");
  params.addCoupledVar("P", -5E7, "The liquid pressure of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("Diffusivity", "The diffusion coefficient of aqeous species");
  
//  params.addRequiredParam<MaterialPropertyName>("D_aq", "The diffusion coefficient of aqueous species");
//  params.addRequiredParam<MaterialPropertyName>("EA_aq", "The activation energy of aqueous species");

  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The alpha value of Van Genuchten equation");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The lambda value of Van Genuchten equation");

  params.addRequiredParam<MaterialPropertyName>("Porosity", "The porosity of porous medium");
  params.addRequiredParam<MaterialPropertyName>("Tortuosity", "The tortuosity of porous medium");

  params.addClassDescription("The diffusion kernel of O2 gaseous / aqueous phase");
  return params;
}

AqPhaseSatDiffusion::AqPhaseSatDiffusion(const InputParameters & parameters)
 : Kernel(parameters),
   _porosity(getMaterialProperty<Real>("Porosity")),
   _tortuosity(getMaterialProperty<Real>("Tortuosity")),
   _T(coupledValue("T")),
   _T_id(coupled("T")),
   _P(coupledValue("P")),
   _P_id(coupled("P")),
   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
//   _D_aq(getMaterialProperty<Real>("D_aq")),
//   _EA_aq(getMaterialProperty<Real>("EA_aq"))
   _Diffusivity(getMaterialProperty<Real>("Diffusivity"))
{
}

Real
AqPhaseSatDiffusion::computeQpResidual()
{
  Real S = pow((1 + pow((-_alpha[_qp] * _P[_qp]), (1 / (1 - _lambda[_qp])))), -_lambda[_qp]);
  Real R = 8.314;
  return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * S * _Diffusivity[_qp] *  _grad_u[_qp];
//  return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * S * _D_aq[_qp] * exp(-_EA_aq[_qp] / (R * _T[_qp])) *  _grad_u[_qp];
}

Real
AqPhaseSatDiffusion::computeQpJacobian()
{
  Real S = pow((1 + pow((-_alpha[_qp] * _P[_qp]), (1 / (1 - _lambda[_qp])))), -_lambda[_qp]);
  Real R = 8.314;
  return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * S * _Diffusivity[_qp] * _grad_phi[_j][_qp];
//  return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * S * _D_aq[_qp] * exp(-_EA_aq[_qp] / (R * _T[_qp])) * _grad_phi[_j][_qp];
}

Real
AqPhaseSatDiffusion::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real S = pow((1 + pow((-_alpha[_qp] * _P[_qp]), (1 / (1 - _lambda[_qp])))), -_lambda[_qp]);
  Real dS_dP = -_lambda[_qp] * pow((1 + pow((-_alpha[_qp] * _P[_qp]), (1 / (1 - _lambda[_qp])))), -_lambda[_qp] - 1) * (_alpha[_qp] / (_lambda[_qp] - 1) * pow((-_alpha[_qp] * _P[_qp]), _lambda[_qp] / (1 - _lambda[_qp])));
  Real R = 8.314;

  if (jvar == _P_id)
      return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * dS_dP * _Diffusivity[_qp] *  _grad_u[_qp];
//      return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * dS_dP * _D_aq[_qp] * exp(-_EA_aq[_qp] / (R * _T[_qp])) *  _grad_u[_qp];
  else if (jvar == _T_id)
      return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * S * _Diffusivity[_qp] *  _grad_u[_qp];
//      return _grad_test[_i][_qp] * _porosity[_qp] * _tortuosity[_qp] * S * _D_aq[_qp] * exp(-_EA_aq[_qp] / (R * _T[_qp])) *  _grad_u[_qp];
  else
      return 0.0;
}
