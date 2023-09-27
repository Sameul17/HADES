# Unsaturated Darcy-Richards flow without using an Action

[Mesh]
  file = '2D_validation_test.msh'
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
    block = 'bentonite backfill hostrock'
  [../]
  [./flux_Pwater]
    type = RichardsFlowAdvectiveFlux
    variable = pwater
    Temp = T
    Water_density = Water_density 
    Water_viscosity = Water_viscosity
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm
    block = 'bentonite backfill hostrock'
  [../]
  [./Gas_Diffusion]
    type = RichardsFlowGasDiffusion
    variable = pwater
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Temp = T
    Water_density = Water_density 
    block = 'bentonite backfill hostrock'
  [../]
 
# Heat Transfer
  [./dT_dt]
    type = RichardsFlowEnergyTimeDerivative
    variable = T
    Medium_density = Medium_density
    Medium_specific_heat = Medium_specific_heat
    block = 'bentonite backfill hostrock'
  [../]
  [./Heat_conduction]
    type = RichardsFlowEnergyConduction
    variable = T
    Thermal_conductivity = Thermal_conductivity
    block = 'bentonite backfill hostrock'
  [../]
  [./Heat_convection]
    type = RichardsFlowEnergyConvection
    variable = T
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm
    Pwater = pwater
    Water_density = Water_density
    Water_viscosity = Water_viscosity
    gravity = '0 0 -9.81'
    block = 'bentonite backfill hostrock'
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
    block = 'bentonite backfill'
   [../]
   [./IC_T_canister]
    type = ConstantIC
    variable = T
    value = 298.15
    boundary = 'spent_fuel spent_fuel_top spent_fuel_bottom'
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
    value = -1.2E7 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_P_backfill]
    type = ConstantIC
    variable = pwater
    value = -4.147E6 #[Pa]
    block = 'backfill'
  [../]
#  [./IC_P_bentonite_surface]
#    type = ConstantIC
#    variable = pwater
#    value = -4207287.175 #[Pa]
#    boundary = 'bentonite_surface'
#  [../]
#  [./IC_P_backfill_surface]
#    type = ConstantIC
#    variable = pwater
#    value = -2272727.273 #[Pa]
#    boundary = 'backfill_surface'
#  [../]
  [./IC_P_hostrock]
    type = ConstantIC
    variable = pwater
    value = -1.549E5
    block = 'hostrock'
  [../]
[]

[BCs]
# Pressure
  [./BC_P]
    type = DirichletBC
    variable = pwater
    value = -1.549E5
    boundary = 'top bottom'
  [../]

# Temperature
  [./BC_T_top]
    type = DirichletBC
    variable = T
    value = '296.65'
    boundary = 'top'
  [../]
  [./BC_T_bottom]
    type = DirichletBC
    variable = T
    value = '299.65'
    boundary = 'bottom'
  [../]
 
# Heat Flux
  [./BC_HeatFlux]
    type = FunctionNeumannBC
    variable = T
    boundary = 'spent_fuel spent_fuel_top spent_fuel_bottom'
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
     vals = '26830 0.51 4.83 31536000'
     value = 'P0 * yts * ((t + 30)^-0.758) / (2 * pi * r * h + 2 * pi * r * r)'
  [../]
[]

[Materials]
  [./permeability_bentonite]
    type = PermeabilityProperty
    block = 'bentonite'
    Pwater = pwater
    n = 1.9
  [../]
  [./permeability_backfill]
    type = PermeabilityProperty
    block = 'backfill'
    Pwater = pwater
    n = 1.9
  [../]
  [./permeability_hostrock]
    type = PermeabilityProperty
    block = 'hostrock'
    Pwater = pwater
    n = 3
  [../]

# Thermal properties
  [./Thermal_conductivity_bentonite]
    type = ThermalProperty
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '2.271E7'
    wet_thermal_conductivity = '3.785E7'
    dry_specific_heat = 966
    wet_specific_heat = 966
    block = 'bentonite'
  [../]
  [./Thermal_conductivity_backfill]
    type = ThermalProperty
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '3.44E7'
    wet_thermal_conductivity = '6.78E7'
    dry_specific_heat = 981
    wet_specific_heat = 981
    block = 'backfill'
  [../]
  [./Thermal_conductivity_hostrock]
    type = ThermalProperty
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '9E7'
    wet_thermal_conductivity = '1E8'
    dry_specific_heat = 820
    wet_specific_heat = 820
    block = 'hostrock'
  [../]

# Hydraulic properties
  [./Bentonite_Hydraulic_Properties]
    type = HydraulicProperty
    Pwater = pwater
    block = 'bentonite'
  [../]
  [./Backfill_Hydraulic_Properties]
    type = HydraulicProperty
    Pwater = pwater
    block = 'backfill'
  [../]
  [./Hostrock_Hydraulic_Properties]
    type = HydraulicProperty
    Pwater = pwater
    block = 'hostrock'
  [../]

  [./Material_properties_bentonite]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm   porosity Medium_density Medium_specific_heat'
    prop_values = '2.6E-7             0.2941               1.5E-21 0.41     2740           966'
    block = 'bentonite'
  [../]
  [./Material_properties_backfill]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm   porosity Medium_density Medium_specific_heat'
    prop_values = '3.3E-7             0.5                  1.6E-20 0.4      2680           981'
    block = 'backfill'
  [../]
  [./Material_properties_hostrock]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm   porosity Medium_density Medium_specific_heat'
    prop_values = '5E-7               0.6                  1.0E-19 0.01     2650           820'
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
  nl_rel_tol = 1E-1

  automatic_scaling = true
  compute_scaling_once = false

  scaling_group_variables = 'T pwater'
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.8
    dt = 0.01
    growth_factor = 1.2
  [../]
[]

[Outputs]
  exodus = true
[]
