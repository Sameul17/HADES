//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PermeabilityPropertyRock.h"

registerMooseObject("HadesApp", PermeabilityPropertyRock);

InputParameters
PermeabilityPropertyRock::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("Pcap","Put the water pressure variable");
  
  params.addParam<Real>("n", 3, "Put the Corey exponent of the phase by using the simple Corey model");
  params.addClassDescription("This material calculates the thermal conductivity along the saturation by using the linear regression");
  return params;
}

PermeabilityPropertyRock::PermeabilityPropertyRock(const InputParameters & parameters)
  : Material(parameters),
    _P(coupledValue("Pcap")),

    _n(getParam<Real>("n")),

    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),

    _Rperm_l(declareProperty<Real>("Rperm_l"))
{
}

void
PermeabilityPropertyRock::initQpStatefulProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  _Rperm_l[_qp] = pow(Sat,0.5) * pow(1 - pow(1 - pow(Sat,1/0.6),0.6),2);
}

void
PermeabilityPropertyRock::computeQpProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  _Rperm_l[_qp] = pow(Sat,0.5) * pow(1 - pow(1 - pow(Sat,1/0.6),0.6),2);
}
