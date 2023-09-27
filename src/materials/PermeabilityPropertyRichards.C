//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PermeabilityPropertyRichards.h"

registerMooseObject("HadesApp", PermeabilityPropertyRichards);

InputParameters
PermeabilityPropertyRichards::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("Pliq","Put the water pressure variable");
  
  params.addParam<Real>("n", 3, "Put the Corey exponent of the phase by using the simple Corey model");
  params.addClassDescription("This material calculates the thermal conductivity along the saturation by using the linear regression");
  return params;
}

PermeabilityPropertyRichards::PermeabilityPropertyRichards(const InputParameters & parameters)
  : Material(parameters),
    _Pliq(coupledValue("Pliq")),

    _n(getParam<Real>("n")),

    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),

    _Rperm_l(declareProperty<Real>("Rperm_l")),
    _Rperm_g(declareProperty<Real>("Rperm_g"))
{
}

void
PermeabilityPropertyRichards::initQpStatefulProperties()
{
  Real Sat = pow((1 + pow((_alpha[_qp] * _Pliq[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  _Rperm_l[_qp] = pow(Sat , _n);
  _Rperm_g[_qp] = 1 - pow(Sat , _n);
}

void
PermeabilityPropertyRichards::computeQpProperties()
{
  Real Sat = pow((1 + pow((_alpha[_qp] * _Pliq[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  _Rperm_l[_qp] = pow(Sat , _n);
  _Rperm_g[_qp] = 1 - pow(Sat , _n);
}
