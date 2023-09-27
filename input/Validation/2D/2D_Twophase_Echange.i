# Unsaturated Darcy-Richards flow without using an Action

[Mesh]
  file = '2D_small_heater_very_course.msh'
  construct_side_list_from_node_list = true
[]

[Variables]
  [./Pcap]
    order = FIRST
    family = LAGRANGE
  [../]
  [./Pgas]
    order = FIRST
    family = LAGRANGE
  [../]
  [./T]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
# Liquid Pressure
  [./dPcap_dt]
    type = TwoPhaseFlowWaterTimeDerivative
    variable = Pcap
    Temp = T
    Pgas = Pgas
    porosity = porosity 
    Water_density = Water_density 
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite hostrock'
  [../]
  [./Liquid_advec_Pcap]
    type = TwoPhaseFlowLiquidAdvectiveFlux
    variable = Pcap
    Pgas = Pgas
    Water_density = Water_density 
    Water_viscosity = Water_viscosity
#    Water_viscosity = 3.171E-11
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    block = 'bentonite hostrock'
  [../]
  [./Vapor_advec_Pcap]
    type = TwoPhaseFlowVaporAdvectiveFlux
    variable = Pcap
    T = T
    Pgas = Pgas
    Water_density = Water_density 
    Gas_viscosity = 5.70776E-13
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability_gas = Rperm_g
    block = 'bentonite hostrock'
  [../]
  [./Gas_Diffusion]
    type = TwoPhaseFlowVaporDiffusion
    variable = Pcap
    Temp = T
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density 
    block = 'bentonite hostrock'
  [../]

# Gas Pressure
  [./dPgas_dt]
    type = TwoPhaseFlowAirTimeDerivative
    variable = Pgas
    Temp = T
    Pcap = Pcap
    porosity = porosity 
    Water_density = Water_density 
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite hostrock'
  [../]
  [./Air_advec_Pcap]
    type = TwoPhaseFlowAirAdvectiveFlux
    variable = Pgas
    T = T
    Pcap = Pcap
    Water_density = Water_density 
    Gas_viscosity = 5.70776E-13
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability_gas = Rperm_g
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite hostrock'
  [../]
  [./Air_Diffusion]
    type = TwoPhaseFlowAirDiffusion
    variable = Pgas
    Temp = T
    Pcap = Pcap
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density 
    block = 'bentonite hostrock'
  [../]
 
# Heat Transfer
  [./dT_dt]
    type = TwoPhaseFlowEnergyTimeDerivative
    variable = T
    Pcap = Pcap
    Pgas = Pgas
    Medium_density = Medium_density
    Medium_specific_heat = Medium_specific_heat
    Water_density = Water_density
    Water_specific_heat = 4.28E3 
    Porosity = porosity
    Gas_specific_heat = 1.01E3
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite hostrock'
  [../]
  [./Heat_conduction]
    type = TwoPhaseFlowEnergyConduction
    variable = T
    Thermal_conductivity = Thermal_conductivity
    block = 'bentonite hostrock'
  [../]
  [./Heat_convection_Liquid]
    type = TwoPhaseFlowEnergyLiquidConvection
    variable = T
    Pcap = Pcap
    Pgas = Pgas
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    Water_density = Water_density
    Water_specific_heat = 4.28E3
#    Water_viscosity = 3.17E-11
    Water_viscosity = Water_viscosity
    gravity = '0 0 -9.81'
    block = 'bentonite hostrock'
  [../]
  [./Heat_convection_gas]
    type = TwoPhaseFlowEnergyGasConvection
    variable = T
    Pcap = Pcap
    Pgas = Pgas
    Intrinsic_permeability = Iperm 
    Relative_permeability_gas = Rperm_g
    Water_density = Water_density
    Water_specific_heat = 4.28E3
    Gas_specific_heat = 1.01E3
    Gas_viscosity = 5.70776E-13
    gravity = '0 0 -9.81'
    block = 'bentonite hostrock'
  [../]
[]

[AuxVariables]
  [./swater]
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
[]

[ICs]
# Temperature
   [./IC_T]
    type = ConstantIC
    variable = T
    value = 298.15
    block = 'bentonite'
   [../]
   [./IC_T_canister]
    type = ConstantIC
    variable = T
    value = 298.15
    boundary = 'heater'
   [../]
   [./IC_T_hostrock]
    type = FunctionIC
    variable = T
    function = underground_temp
    block = 'hostrock'
   [../]

# Liquid Pressure
  [./IC_Pcap_bentonite]
    type = ConstantIC
    variable = Pcap
    value = -79.071E6 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_Pcap_hostrock]
    type = ConstantIC
    variable = Pcap
    value = -2.8733E5
    block = 'hostrock'
  [../]

# Gas pressure
  [./IC_Pgas_bentonite]
    type = ConstantIC
    variable = Pgas
    value = 1E5 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_P_hostrock]
    type = ConstantIC
    variable = Pgas
    value = 1E5
    block = 'hostrock'
  [../]
[]

[BCs]
# Pressure
  [./BC_Pcap]
    type = DirichletBC
    variable = Pcap
    value = -2.8733E5
    boundary = 'top bottom'
  [../]
  [./BC_Pgas_Top]
    type = DirichletBC
    variable = Pgas
    value = 1E5
    boundary = 'top'
  [../]
  [./BC_Pgas_Bottom]
    type = DirichletBC
    variable = Pgas
    value = 1E5
    boundary = 'bottom'
  [../]

# Temperature
  [./BC_T_top]
    type = DirichletBC
    variable = T
    value = '295.15'
    boundary = 'top'
  [../]
  [./BC_T_bottom]
    type = DirichletBC
    variable = T
    value = '301.15'
    boundary = 'bottom'
  [../]

# Heat Flux
  [./BC_HeatFlux]
    type = FunctionNeumannBC
    variable = T
    boundary = 'heater'
    function = Decay_fn
  [../]
[]

[Functions]
  [./underground_temp]
     type = ParsedFunction
     value = '283.15 - 0.03 * y'
  [../]
  [./Decay_fn]
     type = ParsedFunction
     vars = 'P0 r h yts'
     vals = '1127 0.51 4.83 31536000'
     value = 'P0 * yts * ((t + 30)^-0.758)'
  [../]
[]

[Materials]
  [./permeability_bentonite]
    type = PermeabilityPropertyTwoPhase
    block = 'bentonite'
    Pcap = Pcap
    n = 3
  [../]
  [./permeability_hostrock]
    type = PermeabilityPropertyRockTwoPhase
    block = 'hostrock'
    Pcap = Pcap
  [../]

# Thermal properties
  [./Thermal_conductivity_bentonite]
    type = ThermalProperty
    Temp = T
    Pcap = Pcap
    dry_thermal_conductivity = '1.801E7'
    wet_thermal_conductivity = '3.971E7'
    dry_specific_heat = 801.5
    wet_specific_heat = 801.5
    block = 'bentonite'
  [../]
  [./Thermal_conductivity_hostrock]
    type = ThermalProperty
    Temp = T
    Pcap = Pcap
    dry_thermal_conductivity = '9.461E7'
    wet_thermal_conductivity = '9.461E7'
    dry_specific_heat = 900
    wet_specific_heat = 900
    block = 'hostrock'
  [../]

# Hydraulic properties
  [./Material_properties_bentonite]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm   porosity Medium_density Medium_specific_heat'
    prop_values = '2.857E-8           0.3                  2E-21   0.41     1650           966                 '
    block = 'bentonite'
  [../]
  [./Material_properties_hostrock]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm   porosity Medium_density Medium_specific_heat'
    prop_values = '6.803E-7           0.6                  1.0E-17 0.01     2700           900'
    block = 'hostrock'
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
  [./lu]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  [../]
  [./moose]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
    petsc_options_value = 'hypre boomeramg 500'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 1E3
#  nl_abs_tol = 1E-5
  nl_rel_tol = 1E-3

  automatic_scaling = true
  compute_scaling_once = false

  scaling_group_variables = 'T Pcap Pgas'
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.8
    dt = 0.001
    growth_factor = 1.2
  [../]
[]

[Outputs]
  exodus = true
[]
