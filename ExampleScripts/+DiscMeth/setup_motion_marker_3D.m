%======================================================================
%> @file +DiscMeth/setup_motion_marker_3D.m
%> @brief Function to setup motion problem with angle tracking
%>
%> @author Alexander Weiss adapted from Marlies Nitschke
%> @date July, 2025
%======================================================================

% ======================================================================
%> @brief Function to setup motion problem with marker tracking
%>
%> @param   model          Gait3d: Model used for simulation
%> @param   dataFile       String: Filename containg data struct for tracking
%> @param   initialFile    String: Result used for initial guess
%> @param   resultFile     String: Filename to log files
%> @param   W              Struct: Weights for objective terms
%> @param   nNodesEarlier  Double: Number of samples which we start earlier
%> @param   discMeth       String: Name of discretization method (e.g. BE or RIIa-3)
%> @param   nNodes         Int or String: Number of Nodes (e.g. 50 or '50')
%> @retval  problem        Collocation: Generated standing problem
% ======================================================================
function problem = setup_motion_marker_3D(model,dataFile,initialFile,resultFile,W,nNodesEarlier, discMeth, nNodes)

% Load tracking data
trackingData = TrackingData.loadStruct(dataFile);
indicesStartEnd = trackingData.movementEvents.index(strcmp(trackingData.movementEvents.name, 'R_IC'));
trackingData.trimData(indicesStartEnd(1)-nNodesEarlier, indicesStartEnd(2)-1);
N = nNodes;
trackingData.resampleData(N);
%trackingData.nSamples;

% Create problem
Euler = discMeth;
plotLog = 0;
problem = Collocation_RK(model,N,Euler,resultFile,plotLog);

% Add variables which are optimized
states_min = repmat(model.states.xmin,1,N);
states_max = repmat(model.states.xmax,1,N);
xmin_mit = repmat(model.states.xmin, 1, N);
xmax_mit = repmat(model.states.xmax, 1, N);
idxCPxc = model.extractState('xc');
idxCPzc = model.extractState('zc');
p_global_x = [-5, 5];
p_global_z = [-5, 5];
states_min(idxCPxc, :) = p_global_x(1);
states_max(idxCPxc, :) = p_global_x(2);
states_min(idxCPzc, :) = p_global_z(1);
states_max(idxCPzc, :) = p_global_z(2);
problem.addOptimVar('states',states_min,states_max);
problem.addOptimVar('controls',repmat(model.controls.xmin,1,N), repmat(model.controls.xmax,1,N));

% Butcher Tableau
[A, b, c] = DiscMeth.getButcherTableau(discMeth);
for n_k = 1:size(A,1)
    x_states = ['xk', int2str(n_k), '_states'];
    problem.addOptimVar(x_states, xmin_mit, xmax_mit);
end

h = 1/175; % 175 Hz
targetdur =  h*(N-1);
problem.addOptimVar('dur',targetdur,targetdur);

% Initialize the problem with an old result specified in initialFile
problem.makeinitialguess(initialFile); 

% Add tracking terms
trackingDataMar = trackingData.extractData('marker');
trackingDataGRF = trackingData.extractData('GRF', {'GRF_x_r', 'GRF_y_r', 'GRF_z_r', 'GRF_x_l', 'GRF_y_l', 'GRF_z_l'});
problem.addObjective(@trackMarker,W.trackMarker,trackingDataMar);
problem.addObjective(@trackGRF,W.trackGRF,trackingDataGRF);
% Add effort terms
speedWeighting = 0;
problem.addObjective(@effortTermMuscles,W.effMuscles,'volumeweighted',3,speedWeighting);
problem.addObjective(@effortTermTorques,W.effTorques,2, speedWeighting)
% Add regularization term
problem.addObjective(@regTerm,W.reg)

% Add constraints
problem.addConstraint(@dynamicConstraintsRK,repmat(model.constraints.fmin,size(A,1)+1,N-1),repmat(model.constraints.fmax,size(A,1)+1,N-1), discMeth)
problem.addConstraint(@dynamicsFirstNodeConstraint,model.constraints.fmin,model.constraints.fmax)
end
