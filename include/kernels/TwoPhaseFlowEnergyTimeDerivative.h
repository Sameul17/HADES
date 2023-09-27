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
class TwoPhaseFlowEnergyTimeDerivative : public TimeKernel
{
public:
  static InputParameters validParams();

  TwoPhaseFlowEnergyTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Capillary pressure at the nodes
  const VariableValue & _u_old;
  const VariableValue & _Pcap;
  const VariableValue & _Pcap_old;
  unsigned _Pcap_id;
  const VariableValue & _Pgas;
  const VariableValue & _Pgas_old;
  unsigned _Pgas_id;

  /// Density of porous meidum at the nodes,
  const MaterialProperty<Real> & _Mdensity;
  const MaterialProperty<Real> & _MCp;

  /// Density of porous meidum at the nodes,
  const MaterialProperty<Real> & _Wdensity;
  const MaterialProperty<Real> & _WCp;

  /// Density of water at the nodes,
  const MaterialProperty<Real> & _GCp;

  /// Van Genuchten lambda at the nodes,
  const MaterialProperty<Real> & _lambda;

  /// Van Genuchten alpha at the nodes,
  const MaterialProperty<Real> & _alpha;

  /// Porosity at the nodes,
  const MaterialProperty<Real> & _porosity;

};
