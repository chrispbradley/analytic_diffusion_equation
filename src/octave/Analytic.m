
tStart = 0.0;
tStop = 1.0;
#tStep = 0.05;
tStep = 0.01;
degree = 2;

L = 1.0;
numberOfElements = 80;

materialA = 1.0;
materialB = -0.01;

analyticFunction = 1;
analA = 1.0;
analB = 0.0;
analC = 0.0;

function geometricParameters = initialiseGeometricParameters( L, numberOfElements )

  geometricParameters.L = L;
  geometricParameters.numberOfElements = numberOfElements;
  geometricParameters.numberOfNodes = numberOfElements+1;
  
endfunction

function materialParameters = initialiseMaterialParameters( materialConstants )

  materialParameters.a = materialConstants(1);
  materialParameters.b = materialConstants(2);
  
endfunction

function analyticParameters = initialiseAnalyticParameters( analyticFunction, analyticConstants, geometricParameters, materialParameters )

  analyticParameters.analyticFunction = analyticFunction;
  switch (analyticFunction)
    case 1
      analyticParameters.A = analyticConstants(1);
      analyticParameters.B = analyticConstants(2);
      analyticParameters.C = analyticConstants(3);
      analyticParameters.L = geometricParameters.L;
      analyticParameters.sigma = materialParameters.b;
   otherwise
      printf("ERROR: Invalid anlytic function.\n");
      quit
  endswitch

endfunction

function [ dynamicSolver, t ] = initialiseDynamicSolver( tStart, tStop, tStep, degree )

  dynamicSolver.FIRST_DEGREE = 1;
  dynamicSolver.SECOND_DEGREE = 2;
  dynamicSolver.THIRD_DEGREE = 3;
  dynamicSolver.degree = dynamicSolver.FIRST_DEGREE;
  dynamicSolver.startT = 0.0;
  dynamicSolver.sttopT = 1.0;
  dynamicSolver.deltaT = 0.01;
  dynamicSolver.theta1 = 0.0;
  dynamicSolver.theta2 = 0.0;
  dynamicSolver.theta3 = 0.0;
  dynamicSolver.massMatrixFactor = 0.0;
  dynamicSolver.dampingMatrixFactor = 0.0;
  dynamicSolver.stiffnessMatrixFactor = 0.0;
  dynamicSolver.currentFunctionFactor = 0.0;
  dynamicSolver.previousFunctionFactor = 0.0;
  dynamicSolver.previous2FunctionFactor = 0.0;
  dynamicSolver.previous3FunctionFactor = 0.0;

  dynamicSolver.startT = tStart;
  dynamicSolver.stopT = tStop;
  dynamicSolver.deltaT = tStep;
  switch (degree)
    case 1 #Crank-Nicolson
      dynamicSolver.degree = dynamicSolver.FIRST_DEGREE;
      dynamicSolver.theta1 = 0.5;
      dynamicSolver.stiffnessMatrixFactor = dynamicSolver.theta1*dynamicSolver.deltaT;
      dynamicSolver.dampingMatrixFactor = 1.0;
      dynamicSolver.currentFunctionFactor = dynamicSolver.theta1;
      dynamicSolver.previousFunctionFactor = (1.0 - dynamicSolver.theta1);
    case 2 #Newmark
      dynamicSolver.degree = dynamicSolver.SECOND_DEGREE;
      dynamicSolver.theta1 = 2.0;
      dynamicSolver.theta2 = 1.0;
      dynamicSolver.stiffnessMatrixFactor = (dynamicSolver.theta2*dynamicSolver.deltaT*dynamicSolver.deltaT)/2.0;
      dynamicSolver.dampingMatrixFactor = dynamicSolver.theta1*dynamicSolver.deltaT;
      dynamicSolver.massMatrixFactor = 1.0;
      dynamicSolver.currentFunctionFactor = (dynamicSolver.theta2 + dynamicSolver.theta1)/2.0;
      dynamicSolver.previousFunctionFactor = (1.0 - dynamicSolver.theta2);
      dynamicSolver.previous2FunctionFactor = (dynamicSolver.theta2 - dynamicSolver.theta1)/2.0;
    case 3 #Newmark
      dynamicSolver.degree = dynamicSolver.THIRD_DEGREE;
      dynamicSolver.theta1 = 1.0+0.1;
      dynamicSolver.theta2 = 2.0/3.0-0.1+2.0*0.3025;
      dynamicSolver.theta3 = 6.0*0.3025;
      dynamicSolver.stiffnessMatrixFactor = (dynamicSolver.theta3*dynamicSolver.deltaT*dynamicSolver.deltaT*dynamicSolver.deltaT)/6.0;
      dynamicSolver.dampingMatrixFactor = (dynamicSolver.theta2*dynamicSolver.deltaT*dynamicSolver.deltaT)/2.0;
      dynamicSolver.massMatrixFactor = dynamicSolver.theta1*dynamicSolver.deltaT;
      dynamicSolver.currentFunctionFactor = (dynamicSolver.theta3 + 3.0*dynamicSolver.theta2 + 2.0*dynamicSolver.theta1)/6.0;
      dynamicSolver.previousFunctionFactor = (1.0 - (dynamicSolver.theta3 + 2.0*dynamicSolver.theta2 - dynamicSolver.theta1)/2.0);
      dynamicSolver.previous2FunctionFactor = (dynamicSolver.theta3 + dynamicSolver.theta2 - 2.0*dynamicSolver.theta1)/2.0;
      dynamicSolver.previous3FunctionFactor = (dynamicSolver.theta1 - dynamicSolver.theta3)/6.0;
    otherwise
      printf("ERROR: Invalid degree value.\n");
      quit
  endswitch

  t = dynamicSolver.startT;

endfunction

function x = computeGeometry( geometricParameters )

  x = zeros(geometricParameters.numberOfNodes,1);

  x(1) = 0.0;
  for i = 2:geometricParameters.numberOfNodes-1
    x(i) = (i-1)*geometricParameters.L/geometricParameters.numberOfElements;
  endfor
  x(geometricParameters.numberOfNodes) = geometricParameters.L;

endfunction

function [ analyticU, analyticV, analyticA ] = analytic( geometricParameters, analyticParameters, x, t )
  
  analyticU = zeros(geometricParameters.numberOfNodes,1);
  analyticV = zeros(geometricParameters.numberOfNodes,1);
  analyticA = zeros(geometricParameters.numberOfNodes,1);

  switch (analyticParameters.analyticFunction)
    case 1
      gamma = (4.0*pi*pi*analyticParameters.sigma)/(analyticParameters.L*analyticParameters.L);
      delta = (2.0*pi)/analyticParameters.L;
      for i=1:geometricParameters.numberOfNodes
	analyticU(i) = analyticParameters.A*exp(gamma*t)*cos(delta*x(i)+analyticParameters.B)+analyticParameters.C;
        analyticV(i) = gamma*analyticParameters.A*exp(gamma*t)*cos(delta*x(i)+analyticParameters.B);
        analyticA(i) = gamma*gamma*analyticParameters.A*exp(gamma*t)*cos(delta*x(i)+analyticParameters.B);
      endfor
    otherwise
       printf("ERROR: Invalid analytic function value.\n");
       quit
  endswitch

endfunction

function [ currentU, currentV, currentA ] = initialiseProblem( geometricParameters, analyticParameters, dynamicSolver, x, t )
  
  currentU = zeros(geometricParameters.numberOfNodes,1);
  currentV = zeros(geometricParameters.numberOfNodes,1);
  currentA = zeros(geometricParameters.numberOfNodes,1);
  analyticU = zeros(geometricParameters.numberOfNodes,1);
  analyticV = zeros(geometricParameters.numberOfNodes,1);
  analyticA = zeros(geometricParameters.numberOfNodes,1);
  
  [ analyticU, analyticV, analyticA ] = analytic( geometricParameters, analyticParameters, x, t );
  
  switch (dynamicSolver.degree)
    case dynamicSolver.FIRST_DEGREE
      currentU = analyticU;
    case dynamicSolver.SECOND_DEGREE
      currentU = analyticU;
      currentV = analyticV;
    case dynamicSolver.THIRD_DEGREE
      currentU = analyticU;
      currentV = analyticV;
      currentA = analyticA;
    otherwise
      printf("ERROR: Invalid degree value.\n");
      quit
  endswitch

endfunction

function [ M, C, K ] = computeFEMMatrices( geometricParameters, materialParameters )
 
  #For a 1D domain of length L divided into E linear elements then
  #Compute the matrices for a.du/dt + b.d^2u/dx^2 = 0 => C.du/dt + K.u = 0

  K = zeros(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes);
  C = zeros(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes);
  M = zeros(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes);

  K1 = (materialParameters.b*geometricParameters.numberOfElements)/geometricParameters.L;
  K2 = (-1.0*materialParameters.b*geometricParameters.numberOfElements)/geometricParameters.L;
  C1 = (materialParameters.a*geometricParameters.L)/(6.0*geometricParameters.numberOfElements);
  C2 = (materialParameters.a*geometricParameters.L)/(3.0*geometricParameters.numberOfElements);

  K(1,1) = K2;
  K(1,2) = K1;

  C(1,1) = C2;
  C(1,2) = C1;

  for i = 2:geometricParameters.numberOfNodes-1

    K(i,i-1) = K1;
    K(i,i) = 2.0*K2;
    K(i,i+1) = K1;

    C(i,i-1) = C1;
    C(i,i) = 2.0*C2;
    C(i,i+1) = C1;

  endfor

  K(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes-1) = K1;
  K(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes) = K2;

  C(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes-1) = C1;
  C(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes) = C2;

endfunction

function A = computeSolverMatrix( geometricParameters, dynamicSolver, M, C, K )

  A = zeros(geometricParameters.numberOfNodes,geometricParameters.numberOfNodes);
  
  A = dynamicSolver.massMatrixFactor*M + dynamicSolver.dampingMatrixFactor*C + dynamicSolver.stiffnessMatrixFactor*K

  cond(A)
  
endfunction

function [ previousU, previousV, previousA ] = updatePrevious ( geometricParameters, dynamicSolver, currentU, currentV, currentA )

  previousU = zeros(geometricParameters.numberOfNodes,1);
  previousV = zeros(geometricParameters.numberOfNodes,1);
  previousA = zeros(geometricParameters.numberOfNodes,1);
  
  switch (dynamicSolver.degree)
    case dynamicSolver.FIRST_DEGREE
      previousU = currentU;
    case dynamicSolver.SECOND_DEGREE
      previousU = currentU;
      previousV = currentV;
    case dynamicSolver.THIRD_DEGREE
      previousU = currentU;
      previousV = currentV;
      previousA = currentA;
    otherwise
      printf("ERROR: Invalid degree value.\n");
      quit
  endswitch
  
endfunction

function [ currentU, currentV, currentA, alpha ] = updateBoundaryConditions( geometricParameters, dynamicSolver, analyticU, analyticV, analyticA, previousU, previousV, previousA )

  currentU = zeros(geometricParameters.numberOfNodes,1);
  currentV = zeros(geometricParameters.numberOfNodes,1);
  currentA = zeros(geometricParameters.numberOfNodes,1);
  alpha = zeros(geometricParameters.numberOfNodes,1);

  switch (dynamicSolver.degree)
    case dynamicSolver.FIRST_DEGREE
      currentU(1) = analyticU(1);
      currentU(geometricParameters.numberOfNodes) = analyticU(geometricParameters.numberOfNodes);
      alpha(1) = ( currentU(1) - previousU(1) )/dynamicSolver.deltaT;
      alpha(geometricParameters.numberOfNodes) = ( currentU(geometricParameters.numberOfNodes) - previousU(geometricParameters.numberOfNodes) )/dynamicSolver.deltaT;
    case dynamicSolver.SECOND_DEGREE
      currentU(1) = analyticU(1);
      currentU(geometricParameters.numberOfNodes) = analyticU(geometricParameters.numberOfNodes);
      currentV(1) = analyticV(1);
      currentV(geometricParameters.numberOfNodes) = analyticV(geometricParameters.numberOfNodes);
      alpha(1) = ( currentU(1) - previousU(1) - dynamicSolver.deltaT*previousV(1)) *2.0/(dynamicSolver.deltaT*dynamicSolver.deltaT);
      alpha(geometricParameters.numberOfNodes) = ( currentU(geometricParameters.numberOfNodes) - previousU(geometricParameters.numberOfNodes)  - dynamicSolver.deltaT*previousV(geometricParameters.numberOfNodes) )*2.0/(dynamicSolver.deltaT*dynamicSolver.deltaT);
    case dynamicSolver.THIRD_DEGREE
      currentU(1) = analyticU(1);
      currentU(geometricParameters.numberOfNodes) = analyticU(geometricParameters.numberOfNodes);
      currentV(1) = analyticV(1);
      currentV(geometricParameters.numberOfNodes) = analyticV(geometricParameters.numberOfNodes);
      currentA(1) = analyticA(1);
      currentA(geometricParameters.numberOfNodes) = analyticA(geometricParameters.numberOfNodes);
      alpha(1) = ( currentU(1) - previousU(1) - dynamicSolver.deltaT*previousV(1) - (dynamicSolver.deltaT*dynamicSolver.deltaT)/2.0*previousA(1) )*6.0/(dynamicSolver.deltaT*dynamicSolver.deltaT*dynamicSolver.deltaT);
      alpha(1) = ( currentU(geometricParameters.numberOfNodes) - previousU(geometricParameters.numberOfNodes) - dynamicSolver.deltaT*previousV(geometricParameters.numberOfNodes) - (dynamicSolver.deltaT*dynamicSolver.deltaT)/2.0*previousA(geometricParameters.numberOfNodes) )*6.0/(dynamicSolver.deltaT*dynamicSolver.deltaT*dynamicSolver.deltaT);
   otherwise
      printf("ERROR: Invalid degree value.\n");
      quit
  endswitch

endfunction

function [ meanPredictedU, meanPredictedV, meanPredictedA ] = meanPredictedCalculate ( geometricParameters, dynamicSolver, previousU, previousV, previousA )

  meanPredictedU = zeros(geometricParameters.numberOfNodes,1);
  meanPredictedV = zeros(geometricParameters.numberOfNodes,1);
  meanPredictedA = zeros(geometricParameters.numberOfNodes,1);
  
  switch (dynamicSolver.degree)
    case dynamicSolver.FIRST_DEGREE
      meanPredictedU = previousU;
    case dynamicSolver.SECOND_DEGREE
      dynamicSolver.theta1*dynamicSolver.deltaT
      previousU
      previousV
      meanPredictedU = previousU + dynamicSolver.theta1*dynamicSolver.deltaT*previousV;
      meanPredictedV = previousV;
    case dynamicSolver.THIRD_DEGREE
      meanPredictedU = previousU + dynamicSolver.theta1*dynamicSolver.deltaT*previousV + dynamicSolver.theta2*dynamicSolver.deltaT*dynamicSolver.deltaT*previousA/2.0;
      meanPredictedV = previousV + dynamicSolver.theta1*dynamicSolver.deltaT*previousA;
      meanPredictedA = previousA;
    otherwise
      printf("ERROR: Invalid degree value.\n");
      quit
  endswitch
  
endfunction

function beta = computeSolverRHS( geometricParameters, dynamicSolver, meanPredictedU, meanPredictedV, meanPredictedA, M, C, K )

  beta = zeros(geometricParameters.numberOfNodes,1);
  
  beta = -1.0*(M*meanPredictedA + C*meanPredictedV + K*meanPredictedU);
  
endfunction

function [ reducedA, reducedBeta ] = reduceSystem ( geometricParameters, A, alpha, beta )

   reducedA = zeros(geometricParameters.numberOfNodes-1,geometricParameters.numberOfNodes-1);
   reducedBeta = zeros(geometricParameters.numberOfNodes-1,1);

   reducedA = A(2:geometricParameters.numberOfNodes-1,2:geometricParameters.numberOfNodes-1)
   reducedBeta = beta(2:geometricParameters.numberOfNodes-1,1) - A(2:geometricParameters.numberOfNodes-1,1)*alpha(1) - A(2:geometricParameters.numberOfNodes-1,geometricParameters.numberOfNodes)*alpha(geometricParameters.numberOfNodes);

endfunction

function [ solution, residual ] = solveSystem( geometricParameters, A, alpha, beta )

  solution = zeros(geometricParameters.numberOfNodes,1);
  residual = zeros(geometricParameters.numberOfNodes,1);
  reducedA = zeros(geometricParameters.numberOfNodes-1,geometricParameters.numberOfNodes-1);
  reducedAlpha = zeros(geometricParameters.numberOfNodes-1,1);
  reducedBeta = zeros(geometricParameters.numberOfNodes-1,1);

  [ reducedA, reducedBeta ] = reduceSystem( geometricParameters, A, alpha, beta );
  
  reducedAlpha = reducedA\reducedBeta;

  reducedA
  reducedBeta
  reducedAlpha

  solution(1) = alpha(1);
  solution(2:geometricParameters.numberOfNodes-1) = reducedAlpha;
  solution(geometricParameters.numberOfNodes) = alpha(geometricParameters.numberOfNodes);

  residual = A*alpha - beta
  
endfunction

function [ currentU, currentV, currentA ] = currentCalculate ( geometricParameters, dynamicSolver, previousU, previousV, previousA, alpha )

  currentU = zeros(geometricParameters.numberOfNodes,1);
  currentV = zeros(geometricParameters.numberOfNodes,1);
  currentA = zeros(geometricParameters.numberOfNodes,1);
  
  switch (dynamicSolver.degree)
    case dynamicSolver.FIRST_DEGREE
      currentU = previousU + dynamicSolver.deltaT*alpha;
    case dynamicSolver.SECOND_DEGREE
      currentU = previousU + dynamicSolver.deltaT*previousV + dynamicSolver.deltaT*dynamicSolver.deltaT*alpha/2.0;
      currentV = previousV + dynamicSolver.deltaT*alpha;
    case dynamicSolver.THIRD_DEGREE
      currentU = previousU + dynamicSolver.deltaT*previousV + dynamicSolver.deltaT*dynamicSolver.deltaT*previousA/2.0 + dynamicSolver.deltaT*dynamicSolver.deltaT*dynamicSolver.deltaT*alpha/6.0;
      currentV = previousV + dynamicSolver.deltaT*previousA + dynamicSolver.deltaT*dynamicSolver.deltaT*alpha/2.0;
      currentA = previousA + dynamicSolver.deltaT*alpha;
    otherwise
      printf("ERROR: Invalid degree value.\n");
      quit
  endswitch
  
endfunction

function [ err, percentErr, rmsErr ] = error( z, analyticZ )

  [ numRows, numCols ] = size( z );
  err = zeros(numRows,1);
  percentErr = zeros(numRows,1);
  sum = 0.0;
  printf("  row    value analytic      err    %%err \n");
  for i = 1:numRows
    err(i) = (z(i)-analyticZ(i));
    if( abs(analyticZ(i)) > 0.00000001 )
      percentErr(i) = 100.0*(z(i)-analyticZ(i))/analyticZ(i);
    else
      percentErr(i) = 0.0;
    endif
    sum = sum + err(i)*err(i);
    printf("%5d %8.5f %8.5f %8.5f %7.2f\n",i,z(i),analyticZ(i),err(i),percentErr(i));
  endfor
  rmsErr = sqrt(sum/numRows)

endfunction

geometricParameters = initialiseGeometricParameters( L, numberOfElements );

materialParameters = initialiseMaterialParameters( [ materialA, materialB ] );

analyticParameters = initialiseAnalyticParameters( analyticFunction, [ analA, analB, analC ], geometricParameters, materialParameters );

[ dynamicSolver, t ] = initialiseDynamicSolver( tStart, tStop, tStep, degree );

x = computeGeometry( geometricParameters );

[ currentU, currentV, currentA ] = initialiseProblem( geometricParameters, analyticParameters, dynamicSolver, x, t )

[ M, C, K ] = computeFEMMatrices( geometricParameters, materialParameters )

A = computeSolverMatrix( geometricParameters, dynamicSolver, M, C, K );

while t < tStop

  [ previousU, previousV, previousA ] = updatePrevious( geometricParameters, dynamicSolver, currentU, currentV, currentA )

  t = t + tStep

  [ analyticU, analyticV, analyticA ] = analytic( geometricParameters, analyticParameters, x, t )

  [ currentU, currentV, currentA, alpha ] = updateBoundaryConditions( geometricParameters, dynamicSolver, analyticU, analyticV, analyticA, previousU, previousV, previousA )

  [meanPredictedU, meanPredictedV, meanPredictedA ] = meanPredictedCalculate( geometricParameters, dynamicSolver, previousU, previousV, previousA )

  beta = computeSolverRHS( geometricParameters, dynamicSolver, meanPredictedU, meanPredictedV, meanPredictedA, M, C, K );

  [ alpha, r ] = solveSystem( geometricParameters, A,  alpha, beta );

  [ currentU, currentV, currentA ] = currentCalculate( geometricParameters, dynamicSolver,  previousU, previousV, previousA, alpha )

  [ uError, uPerError, uRMSError ] = error( currentU, analyticU );
  [ vError, vPerError, vRMSError ] = error( currentV, analyticV );
  [ aError, aPerError, aRMSError ] = error( currentA, analyticA );

  uRMSError
  vRMSError
  aRMSError

  plot(x, currentU, x, analyticU)
  #plot(x, uError)

  pause(0.50)

  #quit

endwhile

pause(1.0)
