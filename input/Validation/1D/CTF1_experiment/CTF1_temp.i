# Unsaturated Darcy-Richards flow without using an Action

[Mesh]
  file = 'CTF1.msh'
  construct_side_list_from_node_list = true
[]

[Variables]
  [./Pcap]
    order = FIRST
    family = LAGRANGE
  [../]
#  [./Pgas]
#    order = FIRST
#    family = LAGRANGE
#  [../]
  [./T]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
# Water Mass Balance Equation
  [./dPcap_dt]
    type = TwoPhaseFlowWaterTimeDerivative
    variable = Pcap
    Temp = T
#    Pgas = Pgas
    porosity = porosity 
    Water_density = Water_density 
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'matrix'
  [../]
  [./Liquid_advec_Pcap]
    type = TwoPhaseFlowLiquidAdvectiveFlux
    variable = Pcap
#    Pgas = Pgas
    Water_density = Water_density 
    Water_viscosity = 3.171E-11
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    block = 'matrix'
  [../]
  [./Vapor_advec_Pcap]
    type = TwoPhaseFlowVaporAdvectiveFlux
    variable = Pcap
    Temp = T
#    Pgas = Pgas
    Water_density = Water_density 
    Gas_viscosity = 567.648
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability_gas = Rperm_g
    block = 'matrix'
  [../]
  [./Vapor_Diffusion]
    type = TwoPhaseFlowVaporDiffusion
    variable = Pcap
    Temp = T
#    Pgas = Pgas
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density 
    block = 'matrix'
  [../]

# Air Mass Balance Equation
#  [./dPgas_dt]
#    type = TwoPhaseFlowAirTimeDerivative
#    variable = Pgas
#    Temp = T
#    Pcap = Pcap
#    porosity = porosity
#    Water_density = Water_density
#    Van_genuchten_alpha = Van_genuchten_alpha
#    Van_genuchten_lambda = Van_genuchten_lambda
#    block = 'matrix'
#  [../]
#  [./Air_advec_Pgas]
#    type = TwoPhaseFlowAirAdvectiveFlux
#    variable = Pgas
#    Temp = T
#    Pcap = Pcap
#    Water_density = Water_density 
#    Gas_viscosity = 567.648
#    gravity = '0 0 -9.81'
#    Intrinsic_permeability = Iperm 
#    Relative_permeability_gas = Rperm_g
#    block = 'matrix'
#  [../]
#  [./Air_Diffusion]
#    type = TwoPhaseFlowAirDiffusion
#    variable = Pgas
#    Temp = T
#    Pcap = Pcap
#    porosity = porosity
#    Van_genuchten_alpha = Van_genuchten_alpha
#    Van_genuchten_lambda = Van_genuchten_lambda
#    Water_density = Water_density 
#    block = 'matrix'
#  [../]
 
# Heat Transfer
  [./dT_dt]
    type = TwoPhaseFlowEnergyTimeDerivative
    variable = T
    Pcap = Pcap
#    Pgas = Pgas
    Medium_density = Medium_density
    Medium_specific_heat = Medium_specific_heat
    Water_density = Water_density 
    Water_specific_heat = 1000
    Porosity = porosity
    Gas_specific_heat = 1.01E2
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'matrix'
  [../]
  [./Heat_conduction]
    type = TwoPhaseFlowEnergyConduction
    variable = T
   Thermal_conductivity = Thermal_conductivity
    block = 'matrix'
  [../]
  [./Heat_convection_Liquid]
    type = TwoPhaseFlowEnergyLiquidConvection
    variable = T
    Pcap = Pcap
#    Pgas = Pgas
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    Water_density = Water_density
    Water_specific_heat = 1000
    Water_viscosity = 3.171E-11
    gravity = '0 0 -9.81'
    block = 'matrix'
  [../]
  [./Heat_convection_gas]
    type = TwoPhaseFlowEnergyGasConvection
    variable = T
    Pcap = Pcap
#    Pgas = Pgas
    Intrinsic_permeability = Iperm 
    Relative_permeability_gas = Rperm_g
    Water_density = Water_density
    Water_specific_heat = 1000
    Gas_specific_heat = 1.01E3
    Gas_viscosity = 1.8E-5
    gravity = '0 0 -9.81'
    block = 'matrix'
  [../]
[]

[AuxVariables]
  [./swater]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./Iperm]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./swater]
    type = TwoPhaseSatCalculation
    variable = swater
    Pcap = Pcap
  [../]
  [./Iperm_cal]
    type = IpermCalculation
    variable = Iperm
    Pcap = Pcap
  [../]
[]

[ICs]
# Temperature
   [./IC_T_left]
    type = ConstantIC
    variable = T
    value = 373.15
    boundary = 'left'
   [../]
   [./IC_T_right]
    type = ConstantIC
    variable = T
    value = 303.15
    boundary = 'right'
   [../]
   [./IC_T_matrix]
    type = ConstantIC
    variable = T
    value = 303.15
    block = 'matrix'
   [../]

# Pressure
  [./IC_Pcap_matrix]
    type = ConstantIC
    variable = Pcap
    value = -7.5E7 #[Pa]
    block = 'matrix'
  [../]
#  [./IC_Pgas_matrix]
#    type = ConstantIC
#    variable = Pgas
#    value = 1E5 #[Pa]
#    block = 'matrix'
#  [../]
[]

[BCs]
# Temperature
  [./BC_T_top]
    type = DirichletBC
    variable = T
    value = '373.15'
    boundary = 'left'
  [../]
  [./BC_T_bottom]
    type = DirichletBC
    variable = T
    value = '303.15'
    boundary = 'right'
  [../]
[]

[Materials]
  [./permeability_matrix]
    type = PermeabilityPropertyCTF
    block = 'matrix'
    Pcap = Pcap
    n = 12
    k0 = 1E-13
    beta = 4.2
    n0 = 0.44
  [../]

# Thermal properties
  [./Thermal_conductivity_matrix]
    type = ThermalPropertyTwoPhaseCTF
    Temp = T
    Pcap = Pcap
#    Pgas = Pgas
    dry_thermal_conductivity = '1.577E7'
    wet_thermal_conductivity = '4.037E7'
    dry_specific_heat = 966
    wet_specific_heat = 966
    block = 'matrix'
  [../]

  [./Material_properties_matrix]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda porosity Medium_density Medium_specific_heat n0   k0    beta'
    prop_values = '5.556E-8           0.38                 0.44     1650           966                  0.44 1E-13 4.2'
    block = 'matrix'
  [../]
[]

[Preconditioning]
  active = basic
  [./basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  [../]
  [./preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 0.03836 
#  nl_abs_tol = 1E-6
  nl_rel_tol = 1E-3

  automatic_scaling = true
  compute_scaling_once = false

  scaling_group_variables = 'Pcap T'
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.95
    dt = 1E-7
    growth_factor = 1.05
  [../]
[]

[Outputs]
  exodus = true
[]
