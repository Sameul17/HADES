//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThermalPropertyBENT.h"

registerMooseObject("HadesApp", ThermalPropertyBENT);

InputParameters
ThermalPropertyBENT::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("Pcap", -1E5,"Put the water pressure variable");
  params.addCoupledVar("Temp", 300, "Put the temperature variable");
  
  params.addClassDescription("This material calculates the thermal conductivity along the saturation by using the linear regression");
  return params;
}

ThermalPropertyBENT::ThermalPropertyBENT(const InputParameters & parameters)
  : Material(parameters),
    _P(coupledValue("Pcap")),
    _T(coupledValue("Temp")),

    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
    _Mdensity(getMaterialProperty<Real>("Mdensity")),

    _Mspecific_heat(declareProperty<Real>("Mspecific_heat")),
    _Wdensity(declareProperty<Real>("Water_density")),
    _Wviscosity(declareProperty<Real>("Water_viscosity")),
    _Thermal_conductivity(declareProperty<Real>("Thermal_conductivity"))
{
}

void
ThermalPropertyBENT::initQpStatefulProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real StoY = 3600 * 24 * 365;

  _Mspecific_heat[_qp] = (0.5825 * (_Mdensity[_qp]/1000) * 0.2579 * exp(Sat) - 0.2962) * 1000;
  _Wdensity[_qp] = 999.8 + (958.4 - 999.8) * (_T[_qp] - 273.15) / 100;
  _Wviscosity[_qp] = 0.00002414 * std::pow(10, 247.8 / (_T[_qp] - 140)) / StoY;
  _Thermal_conductivity[_qp] = (0.641 * (_Mdensity[_qp]/1000) + 0.624 * Sat - 0.51) * StoY;
}

void
ThermalPropertyBENT::computeQpProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real StoY = 3600 * 24 * 365;

  _Mspecific_heat[_qp] = (0.5825 * (_Mdensity[_qp]/1000) * 0.2579 * exp(Sat) - 0.2962) * 1000;
  _Wdensity[_qp] = 999.8 + (958.4 - 999.8) * (_T[_qp] - 273.15) / 100;
  _Wviscosity[_qp] = 0.00002414 * std::pow(10, 247.8 / (_T[_qp] - 140)) / StoY;
  _Thermal_conductivity[_qp] = (0.641 * (_Mdensity[_qp]/1000) + 0.624 * Sat - 0.51) * StoY;
}
