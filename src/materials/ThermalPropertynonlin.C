//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThermalPropertynonlin.h"

registerMooseObject("HadesApp", ThermalPropertynonlin);

InputParameters
ThermalPropertynonlin::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("Pwater", -1E5,"Put the water pressure variable");
  params.addCoupledVar("Temp", 300, "Put the temperature variable");
  
  params.addParam<Real>("dry_thermal_conductivity", 3, "Put the dry thermal conductivity");
  params.addParam<Real>("wet_thermal_conductivity", 3, "Put the wet thermal conductivity");
  params.addParam<Real>("dry_specific_heat", 3, "Put the dry specific heat");
  params.addParam<Real>("wet_specific_heat", 3, "Put the dry specific heat");
  params.addClassDescription("This material calculates the thermal conductivity along the saturation by using the linear regression");
  return params;
}

ThermalPropertynonlin::ThermalPropertynonlin(const InputParameters & parameters)
  : Material(parameters),
    _P(coupledValue("Pwater")),
    _T(coupledValue("Temp")),

    _dry_thermal_conductivity(getParam<Real>("dry_thermal_conductivity")),
    _wet_thermal_conductivity(getParam<Real>("wet_thermal_conductivity")),
    _dry_specific_heat(getParam<Real>("dry_specific_heat")),
    _wet_specific_heat(getParam<Real>("wet_specific_heat")),

    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),

    _specific_heat(declareProperty<Real>("specific_heat")),
    _Wdensity(declareProperty<Real>("Water_density")),
    _Wviscosity(declareProperty<Real>("Water_viscosity")),
    _Thermal_conductivity(declareProperty<Real>("Thermal_conductivity"))
{
}

void
ThermalPropertynonlin::initQpStatefulProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real StoY = 3600 * 24 * 365;

  _specific_heat[_qp] = _dry_specific_heat + (_dry_specific_heat - _wet_specific_heat) * Sat;
  _Wdensity[_qp] = 999.8 + (958.4 - 999.8) * (_T[_qp] - 273.15) / 100;
  _Wviscosity[_qp] = 1.0016E-3 / StoY;
//  _Wviscosity[_qp] = 0.0000241 * std::pow(10, 247.8 / (_T[_qp] - 140)) / StoY;
  _Thermal_conductivity[_qp] = 1.28 - 0.71 / (1 + pow(2.7182818, (Sat - 0.65)/0.1));
//  _Thermal_conductivity[_qp] = _dry_thermal_conductivity + (_wet_thermal_conductivity - _dry_thermal_conductivity) * Sat;
}

void
ThermalPropertynonlin::computeQpProperties()
{
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real StoY = 3600 * 24 * 365;

  _specific_heat[_qp] = _dry_specific_heat + (_dry_specific_heat - _wet_specific_heat) * Sat;
  _Wdensity[_qp] = 999.8 + (958.4 - 999.8) * (_T[_qp] - 273.15) / 100;
//  _Wviscosity[_qp] = 0.00002414 * std::pow(10, 247.8 / (_T[_qp] - 140)) / StoY;
  _Wviscosity[_qp] = 1.0016E-3 / StoY;
  _Thermal_conductivity[_qp] = 1.28 - 0.71 / (1 + pow(2.7182818,(Sat - 0.65)/0.1));
//  _Thermal_conductivity[_qp] = _dry_thermal_conductivity + (_wet_thermal_conductivity - _dry_thermal_conductivity) * Sat;
}
