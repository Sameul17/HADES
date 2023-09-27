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
class PermeabilityPropertyRockTwoPhase : public Material
{
public:
  static InputParameters validParams();

  PermeabilityPropertyRockTwoPhase(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:
  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _alpha;
  Real _n;

  const VariableValue & _P;

  MaterialProperty<Real> & _Rperm_l;
  MaterialProperty<Real> & _Rperm_g;
};
