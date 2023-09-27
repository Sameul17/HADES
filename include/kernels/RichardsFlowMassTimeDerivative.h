//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeDerivative.h"

/**
 * Kernel = (mass_component - mass_component_old)/dt
 * where mass_component =
 * porosity*sum_phases(density_phase*saturation_phase*massfrac_phase^component)
 * It is lumped to the nodes.
 * If multiply_by_density = false, then density_phase is not included in this above sum, so the
 * Kernel is calculating the time-derivative of fluid volume
 */
class RichardsFlowMassTimeDerivative : public TimeKernel
{
public:
  static InputParameters validParams();

  RichardsFlowMassTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Capillary pressure at the nodes
  const VariableValue & _u_old;

  /// Temperature at the nodes
  const VariableValue & _T;
  const VariableValue & _T_old;
  unsigned _T_id;

  /// Gas pressure at each quadrature point
  const VariableValue & _Pgas;

  /// Porosity at the nodes, but it can depend on grad(variables) which are actually evaluated at the qps
  const MaterialProperty<Real> & _porosity;

  /// Water density at the nodes, but it can depend on grad(variables) which are actually evaluated at the qps
  const MaterialProperty<Real> & _Wdensity;

  /// Old value of porosity
///  const MaterialProperty<Real> & _porosity_old;

  /// d(porosity)/d(PorousFlow variable) - these derivatives will be wrt variables at the nodes
///  const MaterialProperty<std::vector<Real>> & _dporosity_dvar;

  /// d(porosity)/d(grad PorousFlow variable) - remember these derivatives will be wrt grad(vars) at qps
///  const MaterialProperty<std::vector<RealGradient>> & _dporosity_dgradvar;

  /// Van Genuchten parameter 
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _lambda;

};
