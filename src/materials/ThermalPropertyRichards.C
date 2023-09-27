//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ThermalPropertyRichards.h"

registerMooseObject("HadesApp", ThermalPropertyRichards);

InputParameters
ThermalPropertyRichards::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("Pliq", 1E6,"Put the water pressure variable");
  params.addCoupledVar("Temp", 300, "Put the temperature variable");
  
  params.addParam<Real>("dry_thermal_conductivity", 3, "Put the dry thermal conductivity");
  params.addParam<Real>("wet_thermal_conductivity", 3, "Put the wet thermal conductivity");
  params.addParam<Real>("dry_specific_heat", 3, "Put the dry specific heat");
  params.addParam<Real>("wet_specific_heat", 3, "Put the dry specific heat");
  params.addClassDescription("This material calculates the thermal conductivity along the saturation by using the linear regression");
  return params;
}

ThermalPropertyRichards::ThermalPropertyRichards(const InputParameters & parameters)
  : Material(parameters),
    _Pliq(coupledValue("Pliq")),
    _T(coupledValue("Temp")),

    _dry_thermal_conductivity(getParam<Real>("dry_thermal_conductivity")),
    _wet_thermal_conductivity(getParam<Real>("wet_thermal_conductivity")),
    _dry_specific_heat(getParam<Real>("dry_specific_heat")),
    _wet_specific_heat(getParam<Real>("wet_specific_heat")),

    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),

    _Wspecific_heat(declareProperty<Real>("Water_specific_heat")),
    _Wdensity(declareProperty<Real>("Water_density")),
    _Wviscosity(declareProperty<Real>("Water_viscosity")),
    _Thermal_conductivity(declareProperty<Real>("Thermal_conductivity"))
{
}

void
ThermalPropertyRichards::initQpStatefulProperties()
{
  Real Sat = pow((1 + pow((_alpha[_qp] * _Pliq[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real StoY = 3600 * 24 * 365;

  _Wspecific_heat[_qp] = _dry_specific_heat + (_dry_specific_heat - _wet_specific_heat) * Sat;
  _Wdensity[_qp] = 999.8 + (958.4 - 999.8) * (_T[_qp] - 273.15) / 100;
  _Wviscosity[_qp] = 0.00002414 * std::pow(10, 247.8 / (_T[_qp] - 140)) / StoY;
  _Thermal_conductivity[_qp] = _dry_thermal_conductivity + (_wet_thermal_conductivity - _dry_thermal_conductivity) * Sat;
}

void
ThermalPropertyRichards::computeQpProperties()
{
  Real Sat = pow((1 + pow((_alpha[_qp] * _Pliq[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real StoY = 3600 * 24 * 365;

  _Wspecific_heat[_qp] = _dry_specific_heat + (_dry_specific_heat - _wet_specific_heat) * Sat;
  _Wdensity[_qp] = 999.8 + (958.4 - 999.8) * (_T[_qp] - 273.15) / 100;
  _Wviscosity[_qp] = 0.00002414 * std::pow(10, 247.8 / (_T[_qp] - 140)) / StoY;
  _Thermal_conductivity[_qp] = _dry_thermal_conductivity + (_wet_thermal_conductivity - _dry_thermal_conductivity) * Sat;
}
