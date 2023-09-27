//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PermeabilityPropertyTwoPhase.h"

registerMooseObject("HadesApp", PermeabilityPropertyTwoPhase);

InputParameters
PermeabilityPropertyTwoPhase::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("Pcap", -1E5, "Put the water pressure variable");
  
  params.addParam<Real>("n", 3, "Put the Corey exponent of the phase by using the simple Corey model");
  params.addClassDescription("This material calculates the thermal conductivity along the saturation by using the linear regression");
  return params;
}

PermeabilityPropertyTwoPhase::PermeabilityPropertyTwoPhase(const InputParameters & parameters)
  : Material(parameters),
    _P(coupledValue("Pcap")),

    _n(getParam<Real>("n")),

    _b(getMaterialProperty<Real>("Klinkenberg_const")),
    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),

    _Iperm_l(getMaterialProperty<Real>("Iperm_l")),

    _Rperm_l(declareProperty<Real>("Rperm_l")),
//    _Iperm_g(declareProperty<Real>("Iperm_g")),
    _Rperm_g(declareProperty<Real>("Rperm_g"))
{
}

void
PermeabilityPropertyTwoPhase::initQpStatefulProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  _Rperm_l[_qp] = pow(Sat , _n);
//  _Iperm_g[_qp] = _Iperm_l[_qp] * (1 - _b[_qp] / _P[_qp]);
  _Rperm_g[_qp] = 1 - pow(Sat, _n);
}

void
PermeabilityPropertyTwoPhase::computeQpProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  _Rperm_l[_qp] = pow(Sat , _n);
//  _Iperm_g[_qp] = _Iperm_l[_qp] * (1 - _b[_qp] / _P[_qp]);
  _Rperm_g[_qp] = 1 - pow(Sat, _n);
}
