# Unsaturated Darcy-Richards flow without using an Action

[Mesh]
  file = '2D_small_heater.msh'
  construct_side_list_from_node_list = true
[]

[Variables]
  [./pwater]
    order = FIRST
    family = LAGRANGE
  [../]
  [./T]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
# Pressure
  [./dPwater_dt]
    type = RichardsFlowMassTimeDerivative
    Temp = T
    variable = pwater
    porosity = porosity 
    Water_density = Water_density 
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite hostrock'
  [../]
  [./flux_Pwater]
    type = RichardsFlowAdvectiveFlux
    variable = pwater
    Temp = T
    Water_density = Water_density 
    Water_viscosity = Water_viscosity
#    Water_viscosity = 3.17E-11
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    block = 'bentonite hostrock'
  [../]
  [./Gas_Diffusion]
    type = RichardsFlowGasDiffusionTHMC
    variable = pwater
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Temp = T
    Water_density = Water_density 
    block = 'bentonite hostrock'
  [../]
 
# Heat Transfer
  [./dT_dt]
    type = RichardsFlowEnergyTimeDerivative
    variable = T
    Pwater = pwater
    Medium_density = Medium_density
    Medium_specific_heat = Medium_specific_heat
    Water_density = Water_density 
    Porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite hostrock'
  [../]
  [./Heat_conduction]
    type = RichardsFlowEnergyConduction
    variable = T
    Thermal_conductivity = Thermal_conductivity
    block = 'bentonite hostrock'
  [../]
  [./Heat_convection]
    type = RichardsFlowEnergyConvection
    variable = T
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    Pwater = pwater
    Water_density = Water_density
#    Water_viscosity = 3.17E-11
    Water_viscosity = Water_viscosity
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
    type = SatCalculation
    variable = swater
    Pwater = pwater
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

# Pressure
  [./IC_P_bentonite]
    type = ConstantIC
    variable = pwater
    value = -79.071E6 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_P_hostrock]
    type = ConstantIC
    variable = pwater
    value = -2.8733E5
    block = 'hostrock'
  [../]
[]

[BCs]
# Pressure
  [./BC_P]
    type = DirichletBC
    variable = pwater
    value = -2.8733E5
    boundary = 'top bottom'
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
 
#  [./BC_T_heater]
#    type = DirichletBC
#    variable = T
#    value = '353.15'
#    boundary = 'heater'
#  [../]

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
    type = PermeabilityProperty
    block = 'bentonite'
    Pwater = pwater
    n = 3
  [../]
  [./permeability_hostrock]
    type = PermeabilityPropertyRock
    block = 'hostrock'
    Pcap = pwater
  [../]

# Thermal properties
  [./Thermal_conductivity_bentonite]
    type = ThermalPropertyTHMC
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '1.801E7'
    wet_thermal_conductivity = '3.971E7'
    dry_specific_heat = 801.5
    wet_specific_heat = 801.5
    block = 'bentonite'
  [../]
  [./Thermal_conductivity_hostrock]
    type = ThermalProperty
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '9.461E7'
    wet_thermal_conductivity = '9.461E7'
    dry_specific_heat = 900
    wet_specific_heat = 900
    block = 'hostrock'
  [../]

  [./Material_properties_bentonite]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm   porosity Medium_density Medium_specific_heat'
    prop_values = '2.857E-8           0.3                  2E-21   0.41     1600           966'
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
  end_time = 1E3
#  nl_abs_tol = 1E-7
  nl_rel_tol = 1E-3

  automatic_scaling = true
  compute_scaling_once = false

  scaling_group_variables = 'T pwater'
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
