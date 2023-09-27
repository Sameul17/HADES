// 2021.01.24 I have to add reaction products terms. ex) HS- -> Cu+

#include "DiffusionPBC.h"

registerMooseObject("HadesApp", DiffusionPBC);

InputParameters
DiffusionPBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredCoupledVar("Temp","The temperature of each quadrature point");
  params.addRequiredCoupledVar("P","The capillary pressure of each quadrature point");
  params.addRequiredCoupledVar("Reactant1","HS- anions");

  params.addParam<Real>("Stoichiometric_coeff", 1, "The stoichiometric coefficient for calculation of corrosion depth");

  params.addRequiredParam<MaterialPropertyName>("D_aq","The diffusion coefficient of aqeous species");
  params.addRequiredParam<MaterialPropertyName>("EA_aq","The diffusion coefficient of aqeous species");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_alpha","The alpha value");
  params.addRequiredParam<MaterialPropertyName>("Van_genuchten_lambda","The lambda value");
  params.addRequiredParam<MaterialPropertyName>("Porosity","The porosity");
  params.addRequiredParam<MaterialPropertyName>("Tortuosity","The tortuosity");
  params.addClassDescription(
      "Computes a boundary residual contribution consistent with the Diffusion Kernel. "
      "Does not impose a boundary condition; instead computes the boundary "
      "contribution corresponding to the current value of grad(u) and accumulates "
      "it in the residual vector.");
  return params;
}

DiffusionPBC::DiffusionPBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
   _u_old(valueOld()),
   _T(coupledValue("Temp")),
   _T_id(coupled("Temp")),
   _P(coupledValue("P")),
   _P_id(coupled("P")),
   _C1(coupledValue("Reactant1")),
   _C1_id(coupled("Reactant1")),
   _grad_C1(coupledGradient("Reactant1")),

   _coef(getParam<Real>("Stoichiometric_coeff")),

   _D_aq(getMaterialProperty<Real>("D_aq")),
   _EA_aq(getMaterialProperty<Real>("EA_aq")),

   _alpha(getMaterialProperty<Real>("Van_genuchten_alpha")),
   _lambda(getMaterialProperty<Real>("Van_genuchten_lambda")),

   _porosity(getMaterialProperty<Real>("Porosity")),
   _tortuosity(getMaterialProperty<Real>("Tortuosity"))

{
}

Real
DiffusionPBC::computeQpResidual()
{
  Real R = 8.314;
  Real M_Cu = 6.35E-2; //[kg/mol]
  Real rho_Cu = 8960; //[kg/m3]
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));

//  return -_test[_i][_qp] * _D_aq[_qp] * exp(-_EA_aq[_qp] / (R * _T[_qp])) * _grad_C1[_qp] * _normals[_qp] * Sat * _porosity[_qp] * _tortuosity[_qp] * _coef * M_Cu / rho_Cu;
  return -_test[_i][_qp] * _grad_C1[_qp] * _normals[_qp] * _coef * M_Cu / rho_Cu;
}


Real
DiffusionPBC::computeQpJacobian()
{
  return 0;
}

Real
DiffusionPBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real R = 8.314;
  Real M_Cu = 6.35E-2; //[kg/mol]
  Real rho_Cu = 8960; // [kg/m3]
  Real Sat = pow((1 + pow((-_alpha[_qp] * _P[_qp]), (1 / (1 - _lambda[_qp])))),(-_lambda[_qp]));
  Real dSat_dP = _alpha[_qp] * _lambda[_qp] / (1 - _lambda[_qp]) * pow((-_alpha[_qp] * _P[_qp]),(_lambda[_qp] / (1 - _lambda[_qp]))) * pow((1 + pow((-_alpha[_qp] * _P[_qp]),(1 / (1 - _lambda[_qp])))),(-_lambda[_qp] - 1)) * _phi[_j][_qp];

//  Real dm = -_test[_i][_qp] * _D_aq[_qp] * exp(-_EA_aq[_qp] / (R * _T[_qp])) * _grad_C1[_qp] * _normals[_qp] * Sat * _porosity[_qp] * _tortuosity[_qp];
  Real dm = -_test[_i][_qp] * _grad_C1[_qp] * _normals[_qp];

//  Real dm_T  = dm * (_EA_aq[_qp] / (R * pow(_T[_qp],2))) * _phi[_j][_qp] * _coef * M_Cu / rho_Cu;
//  Real dm_P  = dm * dSat_dP / Sat * _coef * M_Cu / rho_Cu;
//  Real dm_C1 = -_test[_i][_qp] * _D_aq[_qp] * exp(-_EA_aq[_qp] / (R * _T[_qp])) * _grad_phi[_j][_qp] * _normals[_qp] * Sat * _porosity[_qp] * _tortuosity[_qp] * _coef * M_Cu / rho_Cu;

  Real dm_T  = dm * _coef * M_Cu / rho_Cu;
  Real dm_P  = dm * _coef * M_Cu / rho_Cu;
  Real dm_C1 = -_test[_i][_qp] * _grad_phi[_j][_qp] * _normals[_qp] * _coef * M_Cu / rho_Cu;

  if (jvar == _T_id)
        return dm_T;
  else if (jvar == _P_id)
        return dm_P;
  else if (jvar == _C1_id)
        return dm_C1;
  else
	return 0.0;  
	
}
