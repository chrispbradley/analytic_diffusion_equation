PROGRAM AnalyticDiffusion
  
  USE OpenCMISS
  
  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters
  REAL(OC_DP), PARAMETER :: MATERIAL_A_PARAM = 1.0_OC_DP
  REAL(OC_DP), PARAMETER :: MATERIAL_SIGMA_PARAM = -0.01_OC_DP

  REAL(OC_DP), PARAMETER :: LENGTH = 1.0_OC_DP

  REAL(OC_DP), PARAMETER :: ANALYTIC_A_PARAM = 1.0_OC_DP
  REAL(OC_DP), PARAMETER :: ANALYTIC_B_PARAM = 0.0_OC_DP
  REAL(OC_DP), PARAMETER :: ANALYTIC_C_PARAM = 0.0_OC_DP

  REAL(OC_DP), PARAMETER :: TIME_START = 0.0_OC_DP
  REAL(OC_DP), PARAMETER :: TIME_STOP = 1.0_OC_DP
  
  INTEGER(OC_Intg), PARAMETER :: OUTPUT_FREQUENCY = 1

  LOGICAL, PARAMETER :: DIRECT_LINEAR_SOLVER = .TRUE.
  
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(OC_Intg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=12
  INTEGER(OC_Intg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=13
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=14
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=15

  !Program variables
  INTEGER(OC_Intg) :: numberOfArguments,argumentLength,status,IO
  INTEGER(OC_Intg) :: dynamicDegree,numberOfGlobalXElements,interpolationType,numberOfGaussXi,numberOfDOFs
  REAL(OC_DP) :: ghostIntegralError(2),ghostRMSError,globalRMSError,integralError(2),localRMSError,timeStep
  CHARACTER(LEN=255) :: commandArgument,filename
  LOGICAL :: directoryExists

  !Program types
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_ControlLoopType) :: controlLoop
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) ::  analyticField,dependentField,equationsSetField,geometricField,materialsField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: dynamicSolver,linearSolver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

  !Generic OpenCMISS variables
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex,err

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalXElements
    IF(numberOfGlobalXElements<=0) CALL HandleError("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) interpolationType
    IF((interpolationType<OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION).OR. &
      & (interpolationType>OC_BASIS_CUBIC_HERMITE_INTERPOLATION)) &
      & CALL HandleError("Invalid Interpolation specification.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) timeStep
    IF(ABS(timeStep)<=1.0E-8) CALL HandleError("Invalid time step.")
    CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 4.")
    READ(commandArgument(1:argumentLength),*) dynamicDegree
    IF((dynamicDegree<OC_SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
      & (dynamicDegree>OC_SOLVER_DYNAMIC_THIRD_DEGREE)) CALL HandleError("Invalid dynamic degree.")
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfGlobalXElements=10
    interpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    timeStep=0.05_OC_DP
    !dynamicDegree=OC_SOLVER_DYNAMIC_FIRST_DEGREE
    dynamicDegree=OC_SOLVER_DYNAMIC_SECOND_DEGREE
  ENDIF
  
  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  WRITE(filename,'(A,"_",I0,"_",I0,"_",F4.2,"_",I0)') "Diffusion",numberOfGlobalXElements,interpolationType,timeStep,dynamicDegree
  CALL OC_OutputSetOn(filename,err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)
  CALL OC_Context_RandomSeedsSet(context,9999,err)

  !Get the computational nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system to be 1D
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,1,err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"DiffusionRegion",err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CALL OC_Basis_NumberOfXiSet(basis,1,err)
  !Set the basis xi interpolation and number of Gauss points
  CALL OC_Basis_InterpolationXiSet(basis,interpolationType,err)
  CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,MAX(3,MIN(4,interpolationType+1)),err)
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(basis,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)
  !Define the mesh on the region
  CALL OC_GeneratedMesh_ExtentSet(generatedMesh,LENGTH,err)
  CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,numberOfGlobalXElements,err)
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)
 
  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------

  !Decompose
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region  
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the scaling to use
  CALL OC_Field_ScalingTypeSet(geometricField,OC_FIELD_NO_SCALING,err)
  !Set the domain to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1,err)
  !Finish creating the field
  CALL OC_Field_CreateFinish(geometricField,err)
  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations_set for a Generalised Burgers's equation
  CALL OC_EquationsSet_Initialise(equationsSet,err)
  CALL OC_Field_Initialise(equationsSetField,err)
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & OC_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,OC_EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER, &
    & equationsSetField,equationsSet,err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIALS FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set material field variables
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  !Finish the equations set material field variables
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)
  !Initialise materials field
  !Set A
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,MATERIAL_A_PARAM,err)
  !Set sigma
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,MATERIAL_SIGMA_PARAM,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set analytic field variables
  CALL OC_Field_Initialise(analyticField,err)
  CALL OC_EquationsSet_AnalyticCreateStart(equationsSet,OC_EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1, &
    & ANALYTIC_FIELD_USER_NUMBER,analyticField,err)
  !Finish the equations set analytic field variables
  CALL OC_EquationsSet_AnalyticCreateFinish(equationsSet,err)
  !Initialise analytic field
  !Set A
  CALL OC_Field_ComponentValuesInitialise(analyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,ANALYTIC_A_PARAM,err)
  !Set B
  CALL OC_Field_ComponentValuesInitialise(analyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,ANALYTIC_B_PARAM,err)
  !Set C
  CALL OC_Field_ComponentValuesInitialise(analyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 3,ANALYTIC_C_PARAM,err)
  !Set L
  CALL OC_Field_ComponentValuesInitialise(analyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 4,LENGTH,err)
  !Set sigma
  CALL OC_Field_ComponentValuesInitialise(analyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 5,MATERIAL_SIGMA_PARAM,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type (Sparse/Full)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_FULL_MATRICES,err)
  !CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output (NoOutput/TimingOutput/MatrixOutput/SolverMatrix/ElementMatrixOutput)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_TIMING,err)
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & OC_PROBLEM_DIFFUSION_EQUATION_TYPE,OC_PROBLEM_LINEAR_DIFFUSION_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(problem,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! CONTROL LOOPS
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem control
  CALL OC_ControlLoop_Initialise(controlLoop,err)
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  !Get the control loop
  CALL OC_Problem_ControlLoopGet(problem,OC_CONTROL_LOOP_NODE,controlLoop,err)
  !Set the times
  CALL OC_ControlLoop_TimesSet(controlLoop,TIME_START,TIME_STOP,timeStep,err)
  !Set the output timing
  CALL OC_ControlLoop_TimeOutputSet(controlLoop,1,err)
  CALL OC_ControlLoop_OutputTypeSet(controlLoop,OC_CONTROL_LOOP_PROGRESS_OUTPUT,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL OC_Problem_SolversCreateStart(problem,err)

  !Get the dymamic solver
  CALL OC_Solver_Initialise(dynamicSolver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,dynamicSolver,err)
  !Set the output type
  !CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_NO_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_MONITOR_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_TIMING_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_SOLVER_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(dynamicSolver,OC_SOLVER_MATRIX_OUTPUT,err)
  !Set dynamic degree
  CALL OC_Solver_DynamicDegreeSet(dynamicSolver,dynamicDegree,err)
  !Set dynamic scheme
  SELECT CASE(dynamicDegree)
  CASE(OC_SOLVER_DYNAMIC_FIRST_DEGREE)
    CALL OC_Solver_DynamicSchemeSet(dynamicSolver,OC_SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err)
  CASE(OC_SOLVER_DYNAMIC_SECOND_DEGREE)
    CALL OC_Solver_DynamicSchemeSet(dynamicSolver,OC_SOLVER_DYNAMIC_NEWMARK1_SCHEME,err)
  CASE(OC_SOLVER_DYNAMIC_THIRD_DEGREE)
    CALL OC_Solver_DynamicSchemeSet(dynamicSolver,OC_SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME,err)
  CASE DEFAULT
    CALL HandleError("Invalid dynamic degree")
  END SELECT

  !Get the dynamic linear solver
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Solver_DynamicLinearSolverGet(dynamicSolver,linearSolver,err)
  !Set the output type
  !CALL OC_Solver_OutputTypeSet(linearSolver,OC_SOLVER_NO_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(linearSolver,OC_SOLVER_MONITOR_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(linearSolver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(linearSolver,OC_SOLVER_TIMING_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(linearSolver,OC_SOLVER_SOLVER_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(linearSolver,OC_SOLVER_MATRIX_OUTPUT,err)
  !Set the solver settings
  IF(DIRECT_LINEAR_SOLVER) THEN
    CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
    CALL OC_Solver_LibraryTypeSet(linearSolver,OC_SOLVER_MUMPS_LIBRARY,err)
  ELSE
    CALL OC_Solver_LinearTypeSet(linearSolver,OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
    CALL OC_Solver_LinearIterativeMaximumIterationsSet(linearSolver,10000,err)
    CALL OC_Solver_LinearIterativeGMRESRestartSet(linearSolver,50,err)
  ENDIF
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(problem,err)


  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the problem solver equations
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  !Get the dynamic solver equations
  CALL OC_Solver_Initialise(dynamicSolver,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,dynamicSolver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Solver_SolverEquationsGet(dynamicSolver,solverEquations,err)
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Set up the boundary conditions
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  !Finish the creation of the equations set boundary conditions
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !-----------------------------------------------------------------------------------------------------------
  ! INITIAL CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Set the initial condition to the analytic solution at the start time
  !Set the analtyic time to the start time
  CALL OC_EquationsSet_AnalyticTimeSet(equationsSet,TIME_START,err)
  !Evaluate the analytic solution
  CALL OC_EquationsSet_AnalyticEvaluate(equationsSet,err)
  !Initialise the dependent field to the analytic
  SELECT CASE(dynamicDegree)
  CASE(OC_SOLVER_DYNAMIC_FIRST_DEGREE)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_ANALYTIC_VALUES_SET_TYPE,1, &
      & dependentfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, &
      & err)
  CASE(OC_SOLVER_DYNAMIC_SECOND_DEGREE)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_ANALYTIC_VALUES_SET_TYPE,1, &
      & dependentfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, &
      & err)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_ANALYTIC_VELOCITY_VALUES_SET_TYPE,1, &
      & dependentfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VELOCITY_VALUES_SET_TYPE,1, &
      & err)
  CASE(OC_SOLVER_DYNAMIC_THIRD_DEGREE)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_ANALYTIC_VALUES_SET_TYPE,1, &
      & dependentfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1, &
      & err)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_ANALYTIC_VELOCITY_VALUES_SET_TYPE,1, &
      & dependentfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VELOCITY_VALUES_SET_TYPE,1, &
      & err)
    CALL OC_Field_ParametersToFieldParametersComponentCopy( &
      & dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_ANALYTIC_ACCELERATION_VALUES_SET_TYPE,1, &
      & dependentfield,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_ACCELERATION_VALUES_SET_TYPE,1, &
      & err)
  CASE DEFAULT
    CALL HandleError("Invalid dynamic degree")
  END SELECT
 
  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------

  INQUIRE(file="./results", exist=directoryExists)
  IF(.NOT.directoryExists) CALL EXECUTE_COMMAND_LINE("mkdir ./output")
  
  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL OC_Problem_Solve(problem,err)
  WRITE(*,'(A)') "Problem solved!"

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Set the analytic time to the stop time
  CALL OC_EquationsSet_AnalyticTimeSet(equationsSet,TIME_STOP,err)
  !Evaluate the analytic solution
  CALL OC_EquationsSet_AnalyticEvaluate(equationsSet,err)
  !Output Analytic analysis
  CALL OC_AnalyticAnalysis_Output(dependentField,filename,err)

  CALL OC_Field_NumberOfDOFsGet(dependentField,OC_FIELD_U_VARIABLE_TYPE,numberOfDOFs,err)
  CALL OC_AnalyticAnalysis_RMSErrorGetNode(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,OC_ANALYTIC_ABSOLUTE_ERROR_TYPE, &
    & localRMSError,ghostRMSError,globalRMSError,err)
  CALL OC_AnalyticAnalysis_IntegralAbsoluteErrorGet(dependentField,OC_FIELD_U_VARIABLE_TYPE,1,integralError, &
    & ghostIntegralError,err)

  OPEN(NEWUNIT=IO,FILE="AnalyticResults.txt",STATUS="UNKNOWN",POSITION="APPEND",ACTION="WRITE")
  WRITE(IO,'(I1," ",I2," ",I1," ",I5," ",E13.6," ",E13.6," ",E13.6," ",E13.6)') interpolationType,numberOfGlobalXElements, &
    & dynamicDegree,numberOfDOFs,LENGTH/REAL(numberOfGlobalXElements),timeStep,globalRMSError,integralError(2)
  CLOSE(IO)

  !export fields
  CALL OC_Fields_Initialise(fields,err)
  CALL OC_Fields_Create(region,fields,err)
  CALL OC_Fields_NodesExport(fields,filename,"FORTRAN",err)
  CALL OC_Fields_ElementsExport(fields,filename,"FORTRAN",err)
  CALL OC_Fields_Finalise(fields,err)
  
  !Destroy the context
  CALL OC_Context_Destroy(context,err)
  !Finalise OpenCMISS
  CALL OC_Finalise(err)
  
  WRITE(*,'(A)') "Program successfully completed."

CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

END PROGRAM AnalyticDiffusion
