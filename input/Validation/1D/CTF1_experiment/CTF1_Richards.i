# Unsaturated Darcy-Richards flow without using an Action

[Mesh]
  file = 'CTF1.msh'
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
    block = 'matrix'
  [../]
  [./flux_Pwater]
    type = RichardsFlowAdvectiveFlux
    variable = pwater
    Temp = T
    Water_density = Water_density 
#    Water_viscosity = Water_viscosity
    Water_viscosity = 3.1710E-11
    gravity = '0 0 -9.81'
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    block = 'matrix'
  [../]
  [./Gas_Diffusion]
    type = RichardsFlowGasDiffusion
    variable = pwater
    porosity = porosity
    Van_genuchten_alpha = Van_genuchten_alpha
    Van_genuchten_lambda = Van_genuchten_lambda
    Temp = T
    Water_density = Water_density 
    block = 'matrix'
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
    block = 'matrix'
  [../]
  [./Heat_conduction]
    type = RichardsFlowEnergyConduction
    variable = T
    Thermal_conductivity = Thermal_conductivity
    block = 'matrix'
  [../]
  [./Heat_convection]
    type = RichardsFlowEnergyConvection
    variable = T
    Intrinsic_permeability = Iperm 
    Relative_permeability = Rperm_l
    Pwater = pwater
    Water_density = Water_density
    Water_viscosity = 3.1710E-11
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
    type = SatCalculation
    variable = swater
    Pwater = pwater
  [../]
  [./Iperm_cal]
    type = IpermCalculation
    variable = Iperm
    Pcap = pwater
  [../]
[]

[ICs]
# Temperature
   [./IC_T]
    type = ConstantIC
    variable = T
    value = 393.15
    boundary = 'left'
   [../]
   [./IC_T_canister]
    type = ConstantIC
    variable = T
    value = 303.15
    boundary = 'right'
   [../]
   [./IC_T_hostrock]
    type = ConstantIC
    variable = T
    value = 303.15
    block = 'matrix'
   [../]

# Pressure
  [./IC_P_matrix]
    type = ConstantIC
    variable = pwater
    value = -75E6 #[Pa]
    block = 'matrix'
  [../]
[]

[BCs]
# Temperature
  [./BC_T_top]
    type = DirichletBC
    variable = T
    value = '393.15'
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
    Pwater = pwater
    n = 12
    k0 = 1E-13
    n0 = 0.44
    beta = 4.2
  [../]

# Thermal properties
  [./Thermal_conductivity_matrix]
    type = ThermalProperty
    Temp = T
    Pwater = pwater
    dry_thermal_conductivity = '1.577E7'
    wet_thermal_conductivity = '4.037E7'
    dry_specific_heat = 966
    wet_specific_heat = 966
    block = 'matrix'
  [../]

# Hydraulic properties
  [./Bentonite_Hydraulic_Properties]
    type = HydraulicProperty
    Pwater = pwater
    block = 'matrix'
  [../]

  [./Material_properties_matrix]
    type = GenericConstantMaterial
    prop_names = 'Van_genuchten_alpha Van_genuchten_lambda porosity Medium_density Medium_specific_heat n0   k0    beta'
    prop_values = '5.5556E-8          0.38                 0.44     1650           966                  0.44 1E-13 4.2'
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
#  nl_abs_tol = 1E-5
  nl_rel_tol = 1E-6

  automatic_scaling = true
  compute_scaling_once = false

  scaling_group_variables = 'T pwater'
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.95
    dt = 1E-9
    growth_factor = 1.05
  [../]
[]

[Outputs]
  exodus = true
[]
