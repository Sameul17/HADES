//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PermeabilityPropertyCTF.h"

registerMooseObject("HadesApp", PermeabilityPropertyCTF);

InputParameters
PermeabilityPropertyCTF::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("Pcap", -1E7,"Put the water pressure variable");

  params.addParam<Real>("k0", 1E-13, "Reference intrinsic permeability");
  params.addParam<Real>("n0", 0.44, "Reference porosity");
  params.addParam<Real>("beta", 0, "Material constant relating to permeability");
  params.addParam<Real>("n", 3, "Put the Corey exponent of the phase by using the simple Corey model");
  params.addClassDescription("This material calculates the thermal conductivity along the saturation by using the linear regression");
  return params;
}

PermeabilityPropertyCTF::PermeabilityPropertyCTF(const InputParameters & parameters)
  : Material(parameters),
    _P(coupledValue("Pcap")),

    _k0(getParam<Real>("k0")),
    _n0(getParam<Real>("n0")),
    _beta(getParam<Real>("beta")),
    _n(getParam<Real>("n")),

    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),

    _Iperm(declareProperty<Real>("Iperm")),
    _Rperm_l(declareProperty<Real>("Rperm_l")),
    _Rperm_g(declareProperty<Real>("Rperm_g"))
{
}

void
PermeabilityPropertyCTF::initQpStatefulProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  Real _nM = _n0 * exp(-_beta * Sat); //Macro porosity
  _Iperm[_qp] = _k0 * pow((1 - _n0) / (1 - _nM) , 2) * pow(_nM / _n0 , 3);
  _Rperm_l[_qp] = pow(Sat , _n);
  _Rperm_g[_qp] = 1 - Sat;
}

void
PermeabilityPropertyCTF::computeQpProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

  Real _nM = _n0 * exp(-_beta * Sat); //Macro porosity
  _Iperm[_qp] = _k0 * pow((1 - _n0) / (1 - _nM) , 2) * pow(_nM / _n0 , 3);
  _Rperm_l[_qp] = pow(Sat , _n);
  _Rperm_g[_qp] = 1 - Sat;
}
