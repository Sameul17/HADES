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
class DiffusionProperty : public Material
{
public:
  static InputParameters validParams();

  DiffusionProperty(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:
  const MaterialProperty<Real> & _porosity;
  const MaterialProperty<Real> & _tortuosity;
  const MaterialProperty<Real> & _lambda;
  const MaterialProperty<Real> & _alpha;


  const VariableValue & _T;
  const VariableValue & _P;

  Real _D_O2_gas;
  Real _EA_O2_gas;
  MaterialProperty<Real> & _RD_O2_gas;
  MaterialProperty<Real> & _REA_O2_gas;

  Real _D_O2_aq;
  Real _EA_O2_aq;
  MaterialProperty<Real> & _RD_O2_aq;
  MaterialProperty<Real> & _REA_O2_aq;

  Real _D_CuCl2_aq;
  Real _EA_CuCl2_aq;
  MaterialProperty<Real> & _RD_CuCl2_aq;
  MaterialProperty<Real> & _REA_CuCl2_aq;

  Real _D_Cu2_aq;
  Real _EA_Cu2_aq;
  MaterialProperty<Real> & _RD_Cu2_aq;
  MaterialProperty<Real> & _REA_Cu2_aq;

  Real _D_Cl_aq;
  Real _EA_Cl_aq;
  MaterialProperty<Real> & _RD_Cl_aq;
  MaterialProperty<Real> & _REA_Cl_aq;

  Real _D_Fe2_aq;
  Real _EA_Fe2_aq;
  MaterialProperty<Real> & _RD_Fe2_aq;
  MaterialProperty<Real> & _REA_Fe2_aq;

  Real _D_HS_aq;
  Real _EA_HS_aq;
  MaterialProperty<Real> & _RD_HS_aq;
  MaterialProperty<Real> & _REA_HS_aq;

  Real _D_SO4_2_aq;
  Real _EA_SO4_2_aq;
  MaterialProperty<Real> & _RD_SO4_2_aq;
  MaterialProperty<Real> & _REA_SO4_2_aq;

  Real _D_H2S_gas;
  Real _EA_H2S_gas;
  MaterialProperty<Real> & _RD_H2S_gas;
  MaterialProperty<Real> & _REA_H2S_gas;

  Real _D_DOM_aq;
  Real _EA_DOM_aq;
  MaterialProperty<Real> & _RD_DOM_aq;
  MaterialProperty<Real> & _REA_DOM_aq;

  MaterialProperty<Real> & _Diffusivity_O2_aq;
  MaterialProperty<Real> & _Diffusivity_O2_two;
  MaterialProperty<Real> & _Diffusivity_HS;
};
