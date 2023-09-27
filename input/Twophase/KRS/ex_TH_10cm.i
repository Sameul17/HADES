# Unsaturated Darcy-Richards flow without using an Action

[Mesh]
  file = '10cm.msh'
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
# Water Mass Balance Equation
  [./dPcap_dt]
    type = TwoPhaseFlowWaterTimeDerivative
    variable = Pcap
    Temp = T
    Pgas = Pgas
    porosity = porosity 
    Water_density = Water_density 
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite backfill hostrock'
  [../]
  [./Liquid_Advec_Pcap]
    type = TwoPhaseFlowLiquidAdvectiveFlux
    variable = Pcap
    Pgas = Pgas
    Water_density = Water_density 
#    Water_viscosity = Water_viscosity
    Water_viscosity = 3.171E-11
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm_l 
    Relative_permeability = Rperm_l
    block = 'bentonite backfill hostrock'
  [../]
  [./Vapor_Advec_Pcap]
    type = TwoPhaseFlowVaporAdvectiveFlux
    variable = Pcap
    Pgas = Pgas
    Temp = T
    Water_density = Water_density 
    Gas_viscosity = 5.70776E-13
    gravity = '0 0 -9.81'
    Intrinsic_permeability_gas = Iperm_g
    Relative_permeability_gas = Rperm_g
#    Van_genuchten_alpha = Van_genuchten_alpha
#    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite backfill hostrock'
  [../]
  [./Gas_Diffusion]
    type = TwoPhaseFlowVaporDiffusion
    variable = Pcap
    Temp = T
    Pgas = Pgas
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density 
    block = 'bentonite backfill hostrock'
  [../]

# Air Mass Balance Equation
  [./dPgas_dt]
    type = TwoPhaseFlowAirTimeDerivative
    variable = Pgas
    Temp = T
    Pcap = Pcap
    porosity = porosity 
    Water_density = Water_density 
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite backfill hostrock'
  [../]
  [./Air_Advec_Pgas]
    type = TwoPhaseFlowAirAdvectiveFlux
    variable = Pgas
    Temp = T
    Pcap = Pcap
    Water_density = Water_density 
    Gas_viscosity = 5.70776E-13
    Intrinsic_permeability_gas = Iperm_g
    Relative_permeability_gas = Rperm_g
    gravity = '0 0 -9.81'
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite backfill hostrock'
  [../]
  [./Air_Diffusion]
    type = TwoPhaseFlowAirDiffusion
    variable = Pgas
    Temp = T
#    Pcap = Pcap
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density 
    block = 'bentonite backfill hostrock'
  [../]
 
# Heat Transfer
  [./dT_dt]
    type = TwoPhaseFlowEnergyTimeDerivative
    variable = T
    Pcap = Pcap
    Pgas = Pgas
    Medium_density = Medium_density
    Medium_specific_heat = Medium_specific_heat
    Porosity = porosity
    Water_specific_heat = 4.28E3
    Gas_specific_heat = 1.01E3
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density
    block = 'bentonite backfill hostrock'
  [../]
  [./Heat_conduction]
    type = TwoPhaseFlowEnergyConduction
    variable = T
    Thermal_conductivity = Thermal_conductivity
    block = 'bentonite backfill hostrock'
  [../]
  [./Heat_convection_liquid]
    type = TwoPhaseFlowEnergyLiquidConvection
    variable = T
    Pcap = Pcap
    Pgas = Pgas
    Intrinsic_permeability = Iperm_l
    Relative_permeability = Rperm_l
    Water_density = Water_density
    Water_specific_heat = 4.28E3
#    Water_viscosity = Water_viscosity
    Water_viscosity = 3.17E-11
    gravity = '0 0 -9.81'
   block = 'bentonite backfill hostrock'
  [../]
  [./Heat_convection_gas]
    type = TwoPhaseFlowEnergyGasConvection
    variable = T
    Pcap = Pcap
    Pgas = Pgas
    Intrinsic_permeability = Iperm_g
    Relative_permeability_gas = Rperm_g
    Water_density = Water_density
    Water_specific_heat = 4.28E3
    Gas_specific_heat = 1.01E3
    Gas_viscosity = 5.70776E-13
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
  [./Swater]
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
  [./IC_Pcap_bentonite]
    type = ConstantIC
    variable = Pcap
    value = -1.143E7 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_Pcap_backfill]
    type = ConstantIC
    variable = Pcap
    value = -4.040E6 #[Pa]
    block = 'backfill'
  [../]
  [./IC_Pcap_hostrock]
    type = ConstantIC
    variable = Pcap
    value = -1.138E5
    block = 'hostrock'
  [../]

  [./IC_Pgas_bentonite]
    type = ConstantIC
    variable = Pgas
    value = 1.01325E5 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_Pgas_backfill]
    type = ConstantIC
    variable = Pgas
    value = 1.01325E5 #[Pa]
    block = 'backfill'
  [../]
  [./IC_Pgas_hostrock]
    type = ConstantIC
    variable = Pgas
    value =  1.01325E5 #[Pa]
    block = 'hostrock'
  [../]
[]

[BCs]
# Pressure
  [./BC_P]
    type = DirichletBC
    variable = Pcap
    value = -1.138E5
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
     expression = '283.15 - 0.03 * y'
  [../]
  [./Decay_fn]
     type = ParsedFunction
     symbol_names = 'P0 r h yts'
     symbol_values = '26830 0.51 4.83 31536000'
     expression = 'P0 * yts * ((t + 40)^-0.758) / (2 * pi * r * h + 2 * pi * r * r)'
  [../]
[]

[Materials]
  [./permeability_bentonite]
    type = PermeabilityPropertyTwoPhase
    block = 'bentonite'
    Pcap = Pcap
    n = 1.9
  [../]
  [./permeability_backfill]
    type = PermeabilityPropertyTwoPhase
    block = 'backfill'
    Pcap = Pcap
    n = 1.9
  [../]
  [./permeability_hostrock]
    type = PermeabilityPropertyTwoPhase
    block = 'hostrock'
    Pcap = Pcap
    n = 3
  [../]

# Thermal properties
  [./Thermal_conductivity_bentonite]
    type = ThermalPropertyTwoPhase
    Temp = T
    Pcap = Pcap
    dry_thermal_conductivity = '2.271E7'
    wet_thermal_conductivity = '3.784E7'
    dry_specific_heat = 966
    wet_specific_heat = 966
    block = 'bentonite'
  [../]
  [./Thermal_conductivity_backfill]
    type = ThermalPropertyTwoPhase
    Temp = T
    Pcap = Pcap
    dry_thermal_conductivity = '3.437E7'
    wet_thermal_conductivity = '6.777E7'
    dry_specific_heat = 981
    wet_specific_heat = 981
    block = 'backfill'
  [../]
  [./Thermal_conductivity_hostrock]
    type = ThermalPropertyTwoPhase
    Temp = T
    Pcap = Pcap
    dry_thermal_conductivity = '8.997E7'
    wet_thermal_conductivity = '9.981E7'
    dry_specific_heat = 820
    wet_specific_heat = 820
    block = 'hostrock'
  [../]

  [./Material_properties_bentonite]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm_l Iperm_g porosity Medium_density Medium_specific_heat Klinkenberg_const'
    prop_values = '2.6E-7             0.2941               1.5E-20 1.0E-14 0.41     2740           966                  1.0E9'
    block = 'bentonite'
  [../]
  [./Material_properties_backfill]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm_l Iperm_g porosity Medium_density Medium_specific_heat Klinkenberg_const'
    prop_values = '3.3E-7             0.5                  1.6E-19 1.0E-13 0.4      2680           981                  1.0E8'
    block = 'backfill'
  [../]
  [./Material_properties_hostrock]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm_l Iperm_g porosity Medium_density Medium_specific_heat Klinkenberg_const'
    prop_values = '5E-7               0.6                  1.0E-18 1.0E-12 0.001    2650           820                  6.86E5'
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
 # nl_abs_tol = 1E-3
  nl_rel_tol = 1E-3

  automatic_scaling = true
  compute_scaling_once = false

  scaling_group_variables = 'Pcap Pgas T'
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 0.001
    growth_factor = 1.11
  [../]
[]

[Outputs]
  exodus = true
[]
