//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"

/**
 * A FluxBC which is consistent with the boundary terms arising from
 * the Diffusion Kernel. The flux vector in this case is simply
 * grad(u) and the residual contribution is:
 *
 * \f$ F(u) = - \int_{\Gamma} \nabla u * \hat n * \phi d\Gamma \f$
 *
 * In contrast to e.g. VectorNeumannBC, the user does not provide a
 * specified value of the flux when using this class, instead the
 * residual contribution corresponding to the current value of grad(u)
 * is computed and accumulated into the residual vector.
 */
class DiffusionPBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  DiffusionPBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _u_old;

  const VariableValue & _T;
  unsigned _T_id;

  const VariableValue & _P;
  unsigned _P_id;

  const VariableValue & _C1;
  const VariableGradient & _grad_C1;
  unsigned _C1_id;

  Real _coef;

  const MaterialProperty<Real> & _D_aq;
  const MaterialProperty<Real> & _EA_aq;

  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> & _lambda;

  const MaterialProperty<Real> & _porosity;
  const MaterialProperty<Real> & _tortuosity;
};
