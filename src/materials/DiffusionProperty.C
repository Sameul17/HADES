//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusionProperty.h"

registerMooseObject("HadesApp", DiffusionProperty);

InputParameters
DiffusionProperty::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("T", 298.15,"Put the temperature variable");
  params.addCoupledVar("P", 1,"Put the pressure variable");

  params.addParam<Real>("D_O2_gas",520,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_O2_gas",-2060,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_O2_aq",5.36E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_O2_aq",15000,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_CuCl2_aq",1.89E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_CuCl2_aq",18800,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_Cu2_aq",1.89E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_Cu2_aq",15000,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_Cl_aq",6.31E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_Cl_aq",15000,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_Fe2_aq",1.58E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_Fe2_aq",15000,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_HS_aq",1.58E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_HS_aq",15000,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_SO4_2_aq",3.15E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_SO4_2_aq",15000,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_H2S_gas",520,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_H2S_gas",-2060,"Put the activation energy of gaseous O2");

  params.addParam<Real>("D_DOM_aq",3.15E-2,"Put the diffusion coefficient of gaseous O2");
  params.addParam<Real>("EA_DOM_aq",15000,"Put the activation energy of gaseous O2");
  return params;
}

DiffusionProperty::DiffusionProperty(const InputParameters & parameters)
  : Material(parameters),
    _T(coupledValue("T")),
    _P(coupledValue("P")),

    _D_O2_gas(getParam<Real>("D_O2_gas")),
    _EA_O2_gas(getParam<Real>("EA_O2_gas")),
    _RD_O2_gas(declareProperty<Real>("RD_O2_gas")),
    _REA_O2_gas(declareProperty<Real>("REA_O2_gas")),

    _D_O2_aq(getParam<Real>("D_O2_aq")),
    _EA_O2_aq(getParam<Real>("EA_O2_aq")),
    _RD_O2_aq(declareProperty<Real>("RD_O2_aq")),
    _REA_O2_aq(declareProperty<Real>("REA_O2_aq")),

    _D_CuCl2_aq(getParam<Real>("D_CuCl2_aq")),
    _EA_CuCl2_aq(getParam<Real>("EA_CuCl2_aq")),
    _RD_CuCl2_aq(declareProperty<Real>("RD_CuCl2_aq")),
    _REA_CuCl2_aq(declareProperty<Real>("REA_CuCl2_aq")),

    _D_Cu2_aq(getParam<Real>("D_Cu2_aq")),
    _EA_Cu2_aq(getParam<Real>("EA_Cu2_aq")),
    _RD_Cu2_aq(declareProperty<Real>("RD_Cu2_aq")),
    _REA_Cu2_aq(declareProperty<Real>("REA_Cu2_aq")),

    _D_Cl_aq(getParam<Real>("D_Cl_aq")),
    _EA_Cl_aq(getParam<Real>("EA_Cl_aq")),
    _RD_Cl_aq(declareProperty<Real>("RD_Cl_aq")),
    _REA_Cl_aq(declareProperty<Real>("REA_Cl_aq")),

    _D_Fe2_aq(getParam<Real>("D_Fe2_aq")),
    _EA_Fe2_aq(getParam<Real>("EA_Fe2_aq")),
    _RD_Fe2_aq(declareProperty<Real>("RD_Fe2_aq")),
    _REA_Fe2_aq(declareProperty<Real>("REA_Fe2_aq")),

    _D_HS_aq(getParam<Real>("D_HS_aq")),
    _EA_HS_aq(getParam<Real>("EA_HS_aq")),
    _RD_HS_aq(declareProperty<Real>("RD_HS_aq")),
    _REA_HS_aq(declareProperty<Real>("REA_HS_aq")),

    _D_SO4_2_aq(getParam<Real>("D_SO4_2_aq")),
    _EA_SO4_2_aq(getParam<Real>("EA_SO4_2_aq")),
    _RD_SO4_2_aq(declareProperty<Real>("RD_SO4_2_aq")),
    _REA_SO4_2_aq(declareProperty<Real>("REA_SO4_2_aq")),

    _D_H2S_gas(getParam<Real>("D_H2S_gas")),
    _EA_H2S_gas(getParam<Real>("EA_H2S_gas")),
    _RD_H2S_gas(declareProperty<Real>("RD_H2S_gas")),
    _REA_H2S_gas(declareProperty<Real>("REA_H2S_gas")),

    _D_DOM_aq(getParam<Real>("D_DOM_aq")),
    _EA_DOM_aq(getParam<Real>("EA_DOM_aq")),
    _RD_DOM_aq(declareProperty<Real>("RD_DOM_aq")),
    _REA_DOM_aq(declareProperty<Real>("REA_DOM_aq")),

   // Diffusivity of O2 and HS- for postprocessor that calculate side flux integral
    _Diffusivity_O2_aq(declareProperty<Real>("Diffusivity_O2_aq")),
    _Diffusivity_O2_two(declareProperty<Real>("Diffusivity_O2_two")),
    _Diffusivity_HS(declareProperty<Real>("Diffusivity_HS")),

    _porosity(getMaterialProperty<Real>("porosity")),
    _tortuosity(getMaterialProperty<Real>("tortuosity")),
    _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),
    _alpha(getMaterialProperty<Real>("Van_genuchten_alpha"))

{
}

void
DiffusionProperty::initQpStatefulProperties()
{
	Real R = 8.314;

///Caculate Diffusion coefficient of each species
//RD: frequency coefficient
//D : diffusion coefficient of species at 25C
	_RD_O2_gas[_qp] = _D_O2_gas * exp(_EA_O2_gas / (R * 298.15));
	_RD_O2_aq[_qp] = _D_O2_aq * exp(_EA_O2_aq / (R * 298.15));
	_RD_CuCl2_aq[_qp] = _D_CuCl2_aq * exp(_EA_CuCl2_aq / (R * 298.15));
	_RD_Cu2_aq[_qp] = _D_Cu2_aq * exp(_EA_Cu2_aq / (R * 298.15));
	_RD_Cl_aq[_qp] = _D_Cl_aq * exp(_EA_Cl_aq / (R * 298.15));
	_RD_Fe2_aq[_qp] = _D_Fe2_aq * exp(_EA_Fe2_aq / (R * 298.15));
	_RD_HS_aq[_qp] = _D_HS_aq * exp(_EA_HS_aq / (R * 298.15));
	_RD_SO4_2_aq[_qp] = _D_SO4_2_aq * exp(_EA_SO4_2_aq / (R * 298.15));
	_RD_H2S_gas[_qp] = _D_H2S_gas * exp(_EA_H2S_gas / (R * 298.15));
	_RD_DOM_aq[_qp] = _D_DOM_aq * exp(_EA_DOM_aq / (R * 298.15));

	_REA_O2_gas[_qp] = _EA_O2_gas;
	_REA_O2_aq[_qp] = _EA_O2_aq;
	_REA_CuCl2_aq[_qp] = _EA_CuCl2_aq;
	_REA_Cu2_aq[_qp] = _EA_Cu2_aq;
	_REA_Cl_aq[_qp] = _EA_Cl_aq;
	_REA_Fe2_aq[_qp] = _EA_Fe2_aq;
	_REA_HS_aq[_qp] = _EA_HS_aq;
	_REA_SO4_2_aq[_qp] = _EA_SO4_2_aq;
	_REA_H2S_gas[_qp] = _EA_H2S_gas;
	_REA_DOM_aq[_qp] = _EA_DOM_aq;

///	_Diffusivity_O2[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_O2_aq[_qp] * exp(-_EA_O2_aq / (R * _T[_qp])) * std::pow(1 + std::pow(-_P[_qp] / _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]) + _RD_O2_gas[_qp] * exp(-_EA_O2_gas / (R * _T[_qp])) * std::pow(1 - std::pow(1 + std::pow(-_P[_qp] / _alpha[_qp],1 / (1 - _lambda[_qp])),-_lambda[_qp]),3) * 13.68 * exp(-_EA_O2_gas / (R * _T[_qp])));
	_Diffusivity_O2_aq[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_O2_aq[_qp] * exp(-_EA_O2_aq / (R * _T[_qp])) * pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]));
	_Diffusivity_O2_two[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_O2_aq[_qp] * exp(-_EA_O2_aq / (R * _T[_qp])) * pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]) + _RD_O2_gas[_qp] * exp(-_EA_O2_gas / (R * _T[_qp])) * pow(1 - pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]), 3));
	_Diffusivity_HS[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_HS_aq[_qp] * exp(-_EA_HS_aq / (R * _T[_qp])) * pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])) ,-_lambda[_qp]));

}

void
DiffusionProperty::computeQpProperties()
{

        Real R = 8.314;

///Caculate Diffusion coefficient of each species
	_RD_O2_gas[_qp] = _D_O2_gas * exp(_EA_O2_gas / (R * 298.15));
	_RD_O2_aq[_qp] = _D_O2_aq * exp(_EA_O2_aq / (R * 298.15));
	_RD_CuCl2_aq[_qp] = _D_CuCl2_aq * exp(_EA_CuCl2_aq / (R * 298.15));
	_RD_Cu2_aq[_qp] = _D_Cu2_aq * exp(_EA_Cu2_aq / (R * 298.15));
	_RD_Cl_aq[_qp] = _D_Cl_aq * exp(_EA_Cl_aq / (R * 298.15));
	_RD_Fe2_aq[_qp] = _D_Fe2_aq * exp(_EA_Fe2_aq / (R * 298.15));
	_RD_HS_aq[_qp] = _D_HS_aq * exp(_EA_HS_aq / (R * 298.15));
	_RD_SO4_2_aq[_qp] = _D_SO4_2_aq * exp(_EA_SO4_2_aq / (R * 298.15));
	_RD_H2S_gas[_qp] = _D_H2S_gas * exp(_EA_H2S_gas / (R * 298.15));
	_RD_DOM_aq[_qp] = _D_DOM_aq * exp(_EA_DOM_aq / (R * 298.15));

	_REA_O2_gas[_qp] = _EA_O2_gas;
	_REA_O2_aq[_qp] = _EA_O2_aq;
	_REA_CuCl2_aq[_qp] = _EA_CuCl2_aq;
	_REA_Cu2_aq[_qp] = _EA_Cu2_aq;
	_REA_Cl_aq[_qp] = _EA_Cl_aq;
	_REA_Fe2_aq[_qp] = _EA_Fe2_aq;
	_REA_HS_aq[_qp] = _EA_HS_aq;
	_REA_SO4_2_aq[_qp] = _EA_SO4_2_aq;
	_REA_H2S_gas[_qp] = _EA_H2S_gas;
	_REA_DOM_aq[_qp] = _EA_DOM_aq;

///	_Diffusivity_O2[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_O2_aq[_qp] * exp(-_EA_O2_aq / (R * _T[_qp])) * std::pow(1 + std::pow(-_P[_qp] / _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]) + _RD_O2_gas[_qp] * exp(-_EA_O2_gas / (R * _T[_qp])) * std::pow(1 - std::pow(1 + std::pow(-_P[_qp] / _alpha[_qp],1 / (1 - _lambda[_qp])),-_lambda[_qp]),3) * 13.68 * exp(-_EA_O2_gas / (R * _T[_qp])));
	_Diffusivity_O2_aq[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_O2_aq[_qp] * exp(-_EA_O2_aq / (R * _T[_qp])) * pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]) );
	_Diffusivity_O2_two[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_O2_aq[_qp] * exp(-_EA_O2_aq / (R * _T[_qp])) * pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]) + _RD_O2_gas[_qp] * exp(-_EA_O2_gas / (R * _T[_qp])) * pow(1 - pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])), -_lambda[_qp]), 3));
	_Diffusivity_HS[_qp] = _porosity[_qp] * _tortuosity[_qp] * (_RD_HS_aq[_qp] * exp(-_EA_HS_aq / (R * _T[_qp])) * pow(1 + pow(-_P[_qp] * _alpha[_qp],1 / (1 - _lambda[_qp])) ,-_lambda[_qp]));
}
