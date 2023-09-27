# Unsaturated Darcy-Richards flow without using an Action
# 3D BMT for the DECOVALEX - III project
# Ref: 
# #1 Numerical study of the THM effects on the near-field safety of a hypothetical nuclear waste repository - BMT1 of the DECOVALEX III project. Part 1: Conceptualization and characterization of the problems and summary of results
# #2 Numerical study of the THM effects on the near-field safety of a hypothetical nuclear waste repository - BMT1 of the DECOVALEX III project. Part 2: Effects of THM coupling in continuous and homogeneous rocks
# #3 A numerical study of THM effects on the near-field safety of a hypothetical nuclear waste repository - BMT1 of the DECOVALEX III project. Part 3: Effects of THM coupling in sparsely fractured rocks

[Mesh]
  file = '3D_BMT_10cm.msh'
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
#    Water_viscosity = Water_viscosity
    Water_viscosity = 3.17E-11
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
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
    Pwater = pwater
    Medium_density = Medium_density
    Medium_specific_heat = Medium_specific_heat
    Porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density 
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
    Relative_permeability = Rperm_l
    Pwater = pwater
    Water_density = Water_density
#    Water_viscosity = Water_viscosity
    Water_viscosity = 3.17E-11
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
    value = -3.963E6 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_P_backfill]
    type = ConstantIC
    variable = pwater
    value = -4.040E6 #[Pa]
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
    value = -3.909E4
    block = 'hostrock'
  [../]
[]

[BCs]
# Pressure
  [./BC_P]
    type = DirichletBC
    variable = pwater
    value = -3.909E4
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
     symbol_names = 'P0 r h yts'
     symbol_values = '12100 0.41 1.73 31536000'
     expression = 'P0 * yts * ((t + 30)^-1.04) / (2 * pi * r * h + 2 * pi * r * r)'
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
    wet_thermal_conductivity = '3.784E7'
    dry_specific_heat = 966
    wet_specific_heat = 966
    block = 'bentonite'
  [../]
  [./Thermal_conductivity_backfill]
    type = ThermalProperty
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '3.437E7'
    wet_thermal_conductivity = '6.777E7'
    dry_specific_heat = 981
    wet_specific_heat = 981
    block = 'backfill'
  [../]
  [./Thermal_conductivity_hostrock]
    type = ThermalProperty
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '8.997E7'
    wet_thermal_conductivity = '9.981E7'
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
    prop_values = '7.289E-7           0.2985               1.6E-20 0.41     1650           426'
    block = 'bentonite'
  [../]
  [./Material_properties_backfill]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm porosity Medium_density Medium_specific_heat'
    prop_values = '3.3E-7             0.5                  1.9E-19 0.333  1800         833'
    block = 'backfill'
  [../]
  [./Material_properties_hostrock]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm   porosity Medium_density Medium_specific_heat'
    prop_values = '5E-7               0.6                  1.0E-18 0.004    1650           820'
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
    cutback_factor = 0.9
    dt = 0.001
    growth_factor = 1.1
  [../]
[]

[Outputs]
  exodus = true
[]
