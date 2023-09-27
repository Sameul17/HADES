
[Problem] #To make cylinderical shape
  coord_type = RZ
[]

[Mesh]
   file = 'ENG10-Hostrock.msh'
   construct_side_list_from_node_list = true
[]

[Variables]
  [./Pcap]
    order = FIRST
    family = LAGRANGE
  [../]

  [./T]
    order = FIRST
    family = LAGRANGE
  [../]

  [./O2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./HS-]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
# Pressure
  [./dPcap_dt]
    type = RichardsFlowMassTimeDerivative
    variable = Pcap
    Temp = T
    porosity = porosity
    Water_density = Water_density
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
  [./Liquid_Advec_Pcap]
    type = RichardsFlowAdvectiveFlux
    variable = Pcap
    Water_density = Water_density
    Water_viscosity = Water_viscosity
#    Water_viscosity = 3.171E-11
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm_l
    Relative_permeability = Rperm_l
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
  [./Gas_Diffusion]
    type = RichardsFlowGasDiffusionTHMC
    variable = Pcap
    Temp = T
    porosity = porosity
    tortuosity = tortuosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]

# Heat Transfer
  [./dT_dt]
    type = RichardsFlowEnergyTimeDerivative
    variable = T
    Pwater = Pcap
    Medium_density = Mdensity 
    Medium_specific_heat = Mspecific_heat
    Porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Water_density = Water_density
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
  [./Heat_conduction]
    type = RichardsFlowEnergyConduction
    variable = T
    Thermal_conductivity = Thermal_conductivity
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
  [./Heat_convection_liquid]
    type = RichardsFlowEnergyConvection
    variable = T
    Pwater = Pcap
    Intrinsic_permeability = Iperm_l
    Relative_permeability = Rperm_l
    Water_density = Water_density
    Water_viscosity = Water_viscosity
#    Water_viscosity = 3.17E-11
    gravity = '0 0 -9.81'
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]

# Chemical transfer
  [./dO2_dt]
    type = TimeDerivative
    variable = O2
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
  [./Diff_O2]
    type = AqPhaseSatDiffusion
    variable = O2
    P = Pcap
    T = T
    Diffusivity = Diffusivity_O2_aq
    Porosity = porosity
    Tortuosity = tortuosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite backfill hostrock1 hostrock2'
  [../] 

  [./dHS-_dt]
    type = TimeDerivative
    variable = HS-
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
  [./Diff_HS-]
    type = AqPhaseSatDiffusion
    variable = HS-
    P = Pcap
    T = T
    Diffusivity = Diffusivity_HS
    Porosity = porosity
    Tortuosity = tortuosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    block = 'bentonite backfill hostrock1 hostrock2'
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
    type = SatCalculation
    variable = swater
    Pwater = Pcap
  [../] 
[]

[ICs]
# Temperature
   [./IC_T_canister]
    type = ConstantIC
    variable = T
    value = 298.15
    boundary = 'canistertop canisterside canisterbottom'
   [../]
   [./IC_T_EBS]
    type = ConstantIC
    variable = T
    value = 298.15
    block = 'bentonite backfill'
   [../]
   [./IC_T_hostrock1hostrock2]
    type = FunctionIC
    variable = T
    function = underground_temp
    block = 'hostrock1 hostrock2'
   [../]

# Pressure
  [./IC_Pcap_bentonite]
    type = ConstantIC
    variable = Pcap
    value = -3.556E7 #[Pa]
    block = 'bentonite'
  [../]
  [./IC_Pcap_backfill]
    type = ConstantIC
    variable = Pcap
    value = -5.152E6 #[Pa]
    block = 'backfill'
  [../]
  [./IC_Pcap_hostrock]
    type = ConstantIC
    variable = Pcap
    value = -1.138E5
    block = 'hostrock1 hostrock2'
  [../]

# Chemicals
  [./IC_O2_EBS]
    type = ConstantIC
    variable = O2
    value = 8.55
    block = 'bentonite backfill'
  [../]
[]

[Materials]
  [./permeability_bentonite]
    type = PermeabilityProperty
    block = 'bentonite'
    Pwater = Pcap
    n = 3
  [../]
  [./permeability_backfill]
    type = PermeabilityProperty
    block = 'backfill'
    Pwater = Pcap
    n = 1.9
  [../]
  [./permeability_hostrock1hostrock2]
    type = PermeabilityProperty
    block = 'hostrock1 hostrock2'
    Pwater = Pcap
    n = 3
  [../]

# Thermal properties
  [./Thermal_conductivity_bentonite]
    type = ThermalPropertyBENT
    Temp = T
    Pcap = Pcap
    block = 'bentonite'
  [../]
  [./Thermal_conductivity_backfill]
    type = ThermalProperty
    Temp = T
    Pcap = Pcap
    dry_specific_heat = 980
    wet_specific_heat = 980
    dry_thermal_conductivity = 1
    wet_thermal_conductivity = 2
    block = 'backfill'
  [../]
  [./Thermal_conductivity_hostrock1hostrock2]
    type = ThermalProperty
    Temp = T
    dry_specific_heat = 820
    wet_specific_heat = 820
    dry_thermal_conductivity = 3.05
    wet_thermal_conductivity = 3.31
    block = 'hostrock1 hostrock2'
  [../]

  [./Material_properties_bentonite]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm_l Iperm_g porosity Mdensity Klinkenberg_const tortuosity'
    prop_values = '1.935E-7           0.2591               9.91E-20 1.0E-14 0.41    1400     1.0E9             0.40'
    block = 'bentonite'
  [../]
  [./Material_properties_backfill]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm_l Iperm_g porosity Mdensity Klinkenberg_const tortuosity'
    prop_values = '3.3E-7             0.5                  1.6E-19 1.0E-13 0.4      1600     1.0E8             0.67'
    block = 'backfill'
  [../]
  [./Material_properties_hostrock1hostrock2]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda Iperm_l Iperm_g porosity Mdensity Klinkenberg_const tortuosity'
    prop_values = '5E-7               0.6                  1.0E-18 1.0E-12 0.001    2650     6.86E5            0.8'
    block = 'hostrock1 hostrock2'
  [../]

  [./Diffusion_coeff]
    type = DiffusionProperty
    P = Pcap
    T = T
    D_O2_aq = 3.15E-1 #[m2/yr]
    D_HS_aq = 3.15E-1 #[m2/yr]
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
[]

[BCs]
# Pressure
  [./BC_P]
    type = DirichletBC
    variable = Pcap
    value = -1.138E5
    boundary = 'hostrocktop hostrockbottom'
  [../]

# Temperature
  [./BC_T_top]
    type = DirichletBC
    variable = T
    value = '296.65'
    boundary = 'hostrocktop'
  [../]
  [./BC_T_bottom]
    type = DirichletBC
    variable = T
    value = '299.65'
    boundary = 'hostrockbottom'
  [../]

# Heat Flux
  [./BC_HeatFlux]
    type = FunctionNeumannBC
    variable = T
    boundary = 'canistertop canisterside canisterbottom'
    function = Decay_fn
  [../]

# Chemicals
  [./BC_O2_canister]
    type = DirichletBC
    variable = O2
    boundary = 'canistertop canisterside canisterbottom'
    value = 0
  [../]
  [./BC_O2_interface]
    type = DirichletBC
    variable = O2
    boundary = 'backfilltop backfillside backfillbottom bentoniteside bentonitebottom'
    value = 0
  [../]

  [./BC_HS-_interface]
    type = DirichletBC
    variable = HS-
    boundary = 'backfilltop backfillside backfillbottom bentoniteside bentonitebottom'
    value = 0.0907
  [../]
  [./BC_HS-_canister]
    type = DirichletBC
    variable = HS-
    boundary = 'canistertop canisterside canisterbottom'
    value = 0
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


[Preconditioning]
  active = basic
  [./basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix -snes_converged_reason'
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
  end_time = 1E7
  nl_rel_tol = 1E-3
#  nl_abs_tol = 1E-5

  automatic_scaling = true
  compute_scaling_once = false

  scaling_group_variables = 'Pcap T'

  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 0.0001
    growth_factor = 1.15
  [../]
[]

[Postprocessors]
  [./O2_side]
    type = SideDiffusiveFluxIntegral
    boundary = 'canisterside'
    variable = O2
    diffusivity = Diffusivity_O2_aq
  [../]
  [./O2_total_side]
    type = TimeIntegratedPostprocessor
    value = O2_side
  [../]
  [./O2_top]
    type = SideDiffusiveFluxIntegral
    boundary = 'canistertop'
    variable = O2
    diffusivity = Diffusivity_O2_aq
  [../]
  [./O2_total_top]
    type = TimeIntegratedPostprocessor
    value = O2_top
  [../]
  [./O2_bottom]
    type = SideDiffusiveFluxIntegral
    boundary = 'canisterbottom'
    variable = O2
    diffusivity = Diffusivity_O2_aq
  [../]
  [./O2_total_bottom]
    type = TimeIntegratedPostprocessor
    value = O2_bottom
  [../]

  [./HS-_side]
    type = SideDiffusiveFluxIntegral
    boundary = 'canisterside'
    variable = HS-
    diffusivity = Diffusivity_HS
  [../]
  [./HS-_total_side]
    type = TimeIntegratedPostprocessor
    value = HS-_side
  [../]
  [./HS-_top]
    type = SideDiffusiveFluxIntegral
    boundary = 'canistertop'
    variable = HS-
    diffusivity = Diffusivity_HS
  [../]
  [./HS-_total_top]
    type = TimeIntegratedPostprocessor
    value = HS-_top
  [../]
  [./HS-_bottom]
    type = SideDiffusiveFluxIntegral
    boundary = 'canisterbottom'
    variable = HS-
    diffusivity = Diffusivity_HS
  [../]
  [./HS-_total_bottom]
    type = TimeIntegratedPostprocessor
    value = HS-_bottom
  [../]

  [./O2]
    type = ElementIntegralVariablePostprocessor
    variable = O2
    block = 'bentonite backfill hostrock1 hostrock2'
  [../]
[]

[Outputs]
  exodus = true
  [./csv]
    type = CSV
  [../]
[]
