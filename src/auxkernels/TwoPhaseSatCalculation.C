//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TwoPhaseSatCalculation.h"

registerMooseObject("HadesApp", TwoPhaseSatCalculation);

InputParameters
TwoPhaseSatCalculation::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("Pcap", "The water pressure of each quadrature point");
  return params;
}

TwoPhaseSatCalculation::TwoPhaseSatCalculation(const InputParameters & parameters)
  : AuxKernel(parameters),
    _Pcap(coupledValue("Pcap")),
    
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda"))
{
}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real
TwoPhaseSatCalculation::computeValue()
{
  return  pow((1 + pow((-_alpha[_qp] * _Pcap[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
}
