#!/usr/bin/env python

import sys

# Intialise OpenCMISS
from opencmiss.opencmiss import OpenCMISS_Python as oc

#-----------------------------------------------------------------------------------------------------------
# SET PROBLEM PARAMETERS
#-----------------------------------------------------------------------------------------------------------

MATERIAL_A_PARAM = 1.0
MATERIAL_SIGMA_PARAM = -1.0

LENGTH = 1.0

ANALYTIC_A_PARAM = 1.0
ANALYTIC_B_PARAM = 1.0
ANALYTIC_C_PARAM = 0.0

TIME_START = 0.0
TIME_STOP = 1.0
TIME_STEP = 0.05

OUTPUT_FREQUENCY = 1

(LINEAR_LAGRANGE_INTERPOLATION,
 QUADRATIC_LAGRANGE_INTERPOLATION,
 CUBIC_LAGRANGE_INTERPOLATION,
 CUBIC_HERMITE_INTERPOLATION) = range(1,5)

(CONTEXT_USER_NUMBER,
 COORDINATE_SYSTEM_USER_NUMBER,
 REGION_USER_NUMBER,
 BASIS_USER_NUMBER,
 GENERATED_MESH_USER_NUMBER,
 MESH_USER_NUMBER,
 DECOMPOSITION_USER_NUMBER,
 DECOMPOSER_USER_NUMBER,
 GEOMETRIC_FIELD_USER_NUMBER,
 EQUATIONS_SET_FIELD_USER_NUMBER,
 DEPENDENT_FIELD_USER_NUMBER,
 MATERIALS_FIELD_USER_NUMBER,
 ANALYTIC_FIELD_USER_NUMBER,
 EQUATIONS_SET_USER_NUMBER,
 PROBLEM_USER_NUMBER) = range(1,16)

numberOfGlobalXElements = 10

interpolationType = LINEAR_LAGRANGE_INTERPOLATION

# Override with command line arguments if need be
if len(sys.argv) > 1:
    if len(sys.argv) > 3:
        sys.exit('ERROR: too many arguments- currently only accepting up to 2 options: numberOfGlobalXElements interpolationType')
    numberOfGlobalXElements = int(sys.argv[1])
    if len(sys.argv) > 2:
        interpolationType = int(sys.argv[2])

if (numberOfGlobalXElements < 0):
    sys.exit('ERROR: number of global X elements must be greater than 0.')
    
if (interpolationType == LINEAR_LAGRANGE_INTERPOLATION):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE
    numberOfGaussXi = 2
elif (interpolationType == QUADRATIC_LAGRANGE_INTERPOLATION):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE
    numberOfGaussXi = 3
elif (interpolationType == CUBIC_LAGRANGE_INTERPOLATION):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.CUBIC_LAGRANGE
    numberOfGaussXi = 4
elif (interpolationType == CUBIC_HERMITE_INTERPOLATION):
    interpolationTypeXi = oc.BasisInterpolationSpecifications.CUBIC_HERMITE
    numberOfGaussXi = 4
else:
    sys.exit('ERROR: invalid interpolation type.')

filenameFormat = "DiffusionEquation_{XElem:0d}_{Interp:0d}"
filename = filenameFormat.format(XElem=numberOfGlobalXElements,Interp=interpolationType)

#-----------------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#-----------------------------------------------------------------------------------------------------------

# Create a context for the example
context = oc.Context()
context.Create(CONTEXT_USER_NUMBER)

oc.OutputSetOn(filename)

# Get the world region
worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

#oc.DiagnosticsSetOn(oc.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["Diffusion_FiniteElementCalculate"])

# Get the computational nodes information
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

#-----------------------------------------------------------------------------------------------------------
# COORDINATE SYSTEM
#-----------------------------------------------------------------------------------------------------------

coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context)
coordinateSystem.DimensionSet(1)
coordinateSystem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# REGION
#-----------------------------------------------------------------------------------------------------------
region = oc.Region()
region.CreateStart(REGION_USER_NUMBER,worldRegion)
region.LabelSet("DiffusionEquation")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BASIS
#-----------------------------------------------------------------------------------------------------------

basis = oc.Basis()
basis.CreateStart(BASIS_USER_NUMBER,context)
basis.TypeSet(oc.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(1)
basis.InterpolationXiSet([interpolationTypeXi])
basis.QuadratureNumberOfGaussXiSet([numberOfGaussXi])
basis.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# MESH
#-----------------------------------------------------------------------------------------------------------

generatedMesh = oc.GeneratedMesh()
generatedMesh.CreateStart(GENERATED_MESH_USER_NUMBER,region)
generatedMesh.TypeSet(oc.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
generatedMesh.ExtentSet([LENGTH])
generatedMesh.NumberOfElementsSet([numberOfGlobalXElements])
mesh = oc.Mesh()
generatedMesh.CreateFinish(MESH_USER_NUMBER,mesh)

#-----------------------------------------------------------------------------------------------------------
# MESH DECOMPOSITION
#-----------------------------------------------------------------------------------------------------------

decomposition = oc.Decomposition()
decomposition.CreateStart(DECOMPOSITION_USER_NUMBER,mesh)
decomposition.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# DECOMPOSER
#-----------------------------------------------------------------------------------------------------------

decomposer = oc.Decomposer()
decomposer.CreateStart(DECOMPOSER_USER_NUMBER,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# GEOMETRIC FIELD
#-----------------------------------------------------------------------------------------------------------

geometricField = oc.Field()
geometricField.CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region)
geometricField.DecompositionSet(decomposition)
geometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

#-----------------------------------------------------------------------------------------------------------
# EQUATION SETS
#-----------------------------------------------------------------------------------------------------------

# Create standard Laplace equations set
equationsSetField = oc.Field()
equationsSet = oc.EquationsSet()
equationsSetSpecification = [oc.EquationsSetClasses.CLASSICAL_FIELD,
                             oc.EquationsSetTypes.DIFFUSION_EQUATION,
                             oc.EquationsSetSubtypes.GENERALISED_DIFFUSION]
equationsSet.CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,
        equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField)
equationsSet.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# DEPENDENT FIELD
#-----------------------------------------------------------------------------------------------------------

dependentField = oc.Field()
equationsSet.DependentCreateStart(DEPENDENT_FIELD_USER_NUMBER,dependentField)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,0.5)

#-----------------------------------------------------------------------------------------------------------
# MATERIALS FIELD
#-----------------------------------------------------------------------------------------------------------

materialsField = oc.Field()
equationsSet.MaterialsCreateStart(MATERIALS_FIELD_USER_NUMBER,materialsField)
materialsField.ComponentInterpolationSet(oc.FieldVariableTypes.U,1,oc.FieldInterpolationTypes.CONSTANT)
materialsField.ComponentInterpolationSet(oc.FieldVariableTypes.U,2,oc.FieldInterpolationTypes.CONSTANT)
equationsSet.MaterialsCreateFinish()

# Initialise materials field
materialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,MATERIAL_A_PARAM)
materialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,MATERIAL_SIGMA_PARAM)

#-----------------------------------------------------------------------------------------------------------
# ANALYTIC FIELD
#-----------------------------------------------------------------------------------------------------------

analyticField = oc.Field()
equationsSet.AnalyticCreateStart(oc.DiffusionAnalyticFunctionTypes.DIFFUSION_EQUATION_ONE_DIM_1,
                                 ANALYTIC_FIELD_USER_NUMBER,analyticField)
equationsSet.AnalyticCreateFinish()

# Set values for analytic field
                                                                                                     
analyticField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,ANALYTIC_A_PARAM)
analyticField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,ANALYTIC_B_PARAM)
analyticField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,3,ANALYTIC_C_PARAM)
analyticField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,4,LENGTH)
analyticField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,5,-MATERIAL_SIGMA_PARAM)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS
#-----------------------------------------------------------------------------------------------------------

equations = oc.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(oc.EquationsSparsityTypes.SPARSE)
#equations.OutputTypeSet(oc.EquationsOutputTypes.NONE)
#equations.OutputTypeSet(oc.EquationsOutputTypes.MATRIX)
equations.OutputTypeSet(oc.EquationsOutputTypes.ELEMENT_MATRIX)
equationsSet.EquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# PROBLEM
#-----------------------------------------------------------------------------------------------------------

problem = oc.Problem()
problemSpecification = [oc.ProblemClasses.CLASSICAL_FIELD,
                        oc.ProblemTypes.DIFFUSION_EQUATION,
                        oc.ProblemSubtypes.LINEAR_DIFFUSION]
problem.CreateStart(PROBLEM_USER_NUMBER,context,problemSpecification)
problem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
# CONTROL LOOPS
#-----------------------------------------------------------------------------------------------------------

# Create control loops
problem.ControlLoopCreateStart()
controlLoop = oc.ControlLoop()
problem.ControlLoopGet([oc.ControlLoopIdentifiers.NODE],controlLoop)
controlLoop.TimesSet(TIME_START,TIME_STOP,TIME_STEP)
controlLoop.OutputTypeSet(oc.ControlLoopOutputTypes.TIMING)
controlLoop.TimeOutputSet(OUTPUT_FREQUENCY)
problem.ControlLoopCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVER
#-----------------------------------------------------------------------------------------------------------

# Create problem solver
dynamicSolver = oc.Solver()
problem.SolversCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,dynamicSolver)
#dynamicSolver.OutputTypeSet(oc.SolverOutputTypes.SOLVER)
dynamicSolver.OutputTypeSet(oc.SolverOutputTypes.MATRIX)
#dynamicSolver.LinearTypeSet(oc.LinearSolverTypes.ITERATIVE)
#dynamicSolver.LinearIterativeAbsoluteToleranceSet(1.0E-12)
#dynamicSolver.LinearIterativeRelativeToleranceSet(1.0E-12)
problem.SolversCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVER EQUATIONS
#-----------------------------------------------------------------------------------------------------------

# Create solver equations and add equations set to solver equations
dynamicSolver = oc.Solver()
solverEquations = oc.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,dynamicSolver)
dynamicSolver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(oc.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------------------------

# Create analytic boundary conditions
boundaryConditions = oc.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
solverEquations.BoundaryConditionsAnalytic()
solverEquations.BoundaryConditionsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
# SOLVE
#-----------------------------------------------------------------------------------------------------------

problem.Solve()

#-----------------------------------------------------------------------------------------------------------
# ANALYTIC ANALYSIS
#-----------------------------------------------------------------------------------------------------------

equationsSet.AnalyticTimeSet(TIME_STOP)
oc.AnalyticAnalysis_Output(dependentField,filename)

#-----------------------------------------------------------------------------------------------------------
# OUTPUT
#-----------------------------------------------------------------------------------------------------------

# Export results
fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("DiffusionEquation","FORTRAN")
fields.ElementsExport("DiffusionEquation","FORTRAN")
fields.Finalise()

# Finalise OpenCMISS
oc.Finalise()
