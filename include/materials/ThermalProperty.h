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
class ThermalProperty : public Material
{
public:
  static InputParameters validParams();

  ThermalProperty(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:
  const VariableValue & _P;
  const VariableValue & _T;

  Real _Kdry;
  Real _Kwet;
  Real _Cpdry;
  Real _Cpwet;

  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _Mdensity;

  MaterialProperty<Real> & _Mspecific_heat;
  MaterialProperty<Real> & _Wdensity;
  MaterialProperty<Real> & _Wviscosity;
  MaterialProperty<Real> & _Thermal_conductivity;

};
