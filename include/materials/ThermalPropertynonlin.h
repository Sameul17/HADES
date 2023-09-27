//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "PorousFlowDictator.h"

/**
 * Stateful material class that defines a few properties.
 */
class ThermalPropertynonlin : public Material
{
public:
  static InputParameters validParams();

  ThermalPropertynonlin(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:
  const VariableValue & _P;
  const VariableValue & _T;

  Real _dry_thermal_conductivity;
  Real _wet_thermal_conductivity;
  Real _dry_specific_heat;
  Real _wet_specific_heat;

  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _alpha;

  MaterialProperty<Real> & _specific_heat;
  MaterialProperty<Real> & _Wdensity;
  MaterialProperty<Real> & _Wviscosity;
  MaterialProperty<Real> & _Thermal_conductivity;

};
