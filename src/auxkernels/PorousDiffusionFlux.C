//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousDiffusionFlux.h"
#include "Assembly.h"

registerMooseObject("HadesApp", PorousDiffusionFlux);


InputParameters
PorousDiffusionFlux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Compute components of flux vector for diffusion problems "
                             "$(\\vv{J} = -D \\nabla C)$.");
  params.addRequiredCoupledVar("diffusion_variable", "The name of the variable");
  params.addRequiredCoupledVar("Pwater", "The water pressure of each quadrature point");
  params.addRequiredCoupledVar("Temp", "The temperature of each quadrature point");
  params.addRequiredParam<MaterialPropertyName>("D_aq", "The diffusion coefficient of aqeous species");
  params.addRequiredParam<MaterialPropertyName>("EA_aq", "The activation energy of aqeous species");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha", "The activation energy of aqeous species");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda", "The activation energy of aqeous species");
  params.addRequiredParam<MaterialPropertyName>("porosity", "The activation energy of aqeous species");
  params.addRequiredParam<MaterialPropertyName>("tortuosity", "The activation energy of aqeous species");

  return params;
}

PorousDiffusionFlux::PorousDiffusionFlux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad_u(coupledGradient("diffusion_variable")),
    _P(coupledValue("Pwater")),
    _T(coupledValue("Temp")),
    _porosity(getMaterialProperty<Real>("porosity")),
    _tortuosity(getMaterialProperty<Real>("tortuosity")),
    _D_aq(getMaterialProperty<Real>("D_aq")),
    _EA_aq(getMaterialProperty<Real>("EA_aq")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _normals(_assembly.normals())
{
}

Real
PorousDiffusionFlux::computeValue()
{
  Real gradient = _grad_u[_qp] * _normals[_qp];
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real R = 8.314;

  Real Flux = gradient *  _porosity[_qp] * _tortuosity[_qp] * _D_aq[_qp] * exp(-_EA_aq[_qp] / R * _T[_qp]) * Sat ;
  return Flux;
}
