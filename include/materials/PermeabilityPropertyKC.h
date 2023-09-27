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

/**
 * Stateful material class that defines a few properties.
 */
class PermeabilityPropertyKC : public Material
{
public:
  static InputParameters validParams();

  PermeabilityPropertyKC(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:
  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _alpha;
  Real _k0;
  Real _n0;
  Real _beta;
  Real _n;

  const VariableValue & _P;

  MaterialProperty<Real> & _Iperm;
  MaterialProperty<Real> & _Rperm_l;
  MaterialProperty<Real> & _Rperm_g;
};
