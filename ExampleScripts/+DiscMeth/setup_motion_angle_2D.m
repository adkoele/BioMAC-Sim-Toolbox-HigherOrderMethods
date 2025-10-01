%======================================================================
%> @file +DiscMeth/setup_motion_angle_2D.m
%> @brief Function to setup motion problem with angle tracking
%>
%> @author Alexander Weiss adapted from Marlies Nitschke
%> @date July, 2025
%======================================================================

% ======================================================================
%> @brief Function to setup motion problem with angle tracking
%>
%> @param   model          Gait3d: Model used for simulation
%> @param   dataInvFile    String: Filename containg data struct for angles for tracking
%> @param   dataMeasFile   String: Filename containg data struct for GRFs for tracking
%> @param   initialFile    String: Result used for initial guess
%> @param   resultFile     String: Filename to log files
%> @param   W              Struct: Weights for objective terms
%> @param   nNodesEarlier  Double: Number of samples which we start earlier
%> @param   discMeth       String: Name of discretization method (e.g. BE or RIIa-3)
%> @param   nNodes         Int or String: Number of Nodes (e.g. 50 or '50')
%> @retval  problem        Collocation: Generated standing problem
% ======================================================================
function problem = setup_motion_angle_2D(model,dataInvFile,dataMeasFile,initialFile,resultFile,W, nNodesEarlier, discMeth, nNodes)

N = nNodes;
% Load tracking data
trackingDataInv = TrackingData.loadStruct(dataInvFile);
indicesStartEnd = trackingDataInv.movementEvents.index(strcmp(trackingDataInv.movementEvents.name, 'R_IC'));
trackingDataInv.trimData(indicesStartEnd(1)-nNodesEarlier, indicesStartEnd(2)-1);
trackingDataInv.resampleData(N);

trackingDataMeas = TrackingData.loadStruct(dataMeasFile);
indicesStartEnd = trackingDataMeas.movementEvents.index(strcmp(trackingDataMeas.movementEvents.name, 'R_IC'));
trackingDataMeas.trimData(indicesStartEnd(1)-nNodesEarlier, indicesStartEnd(2)-1);
initSamples = trackingDataMeas.nSamples;
trackingDataMeas.resampleData(N);
%trackingDataMeas.nSamples;

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
p_global_x = [-5, 5];
states_min(idxCPxc, :) = p_global_x(1);
states_max(idxCPxc, :) = p_global_x(2);

problem.addOptimVar('states',states_min,states_max);
problem.addOptimVar('controls',repmat(model.controls.xmin,1,N), repmat(model.controls.xmax,1,N));

% Butcher Tableau
[A, b, c] = DiscMeth.getButcherTableau(discMeth);
for n_k = 1:size(A,1)
    x_states = ['xk', int2str(n_k), '_states'];
    problem.addOptimVar(x_states, xmin_mit, xmax_mit);
end

h = 1/175; % 175 Hz
targetdur =  h*(initSamples-1);
problem.addOptimVar('dur',targetdur,targetdur);

% Initialize the problem with an old result specified in initialFile
problem.makeinitialguess(initialFile); 

% Add tracking terms
trackingDataTrans = trackingDataInv.extractData('translation', {'pelvis_tx', 'pelvis_ty'});
trackingDataAng = trackingDataInv.extractData('angle', {'pelvis_tilt' ...
    'hip_flexion_r' 'knee_angle_r'    'ankle_angle_r' ...
    'hip_flexion_l' 'knee_angle_l'    'ankle_angle_l' });
trackingDataGRF = trackingDataMeas.extractData('GRF', {'GRF_x_r', 'GRF_y_r', 'GRF_x_l', 'GRF_y_l'});
%problem.addObjective(@trackTranslations,W.trackTrans,trackingDataTrans);
problem.addObjective(@trackAngles,W.trackAngle,trackingDataAng);
problem.addObjective(@trackGRF,W.trackGRF,trackingDataGRF);

indTransl = find(strcmp(trackingDataTrans.variables.name, 'pelvis_tx'));
targetSpeed = ((trackingDataTrans.variables.mean{indTransl,1}(end)*1e-3 - trackingDataTrans.variables.mean{indTransl,1}(1)*1e-3) / targetdur);
%problem.addOptimVar('speed',targetSpeed, targetSpeed);

% Add effort terms
speedWeighting = 0;
problem.addObjective(@effortTermMuscles,W.effMuscles,'equal',3,speedWeighting);
% add regularization term
problem.addObjective(@regTerm,W.reg)

isSymmetric = 0;
% Add constraints

%problem.addConstraint(@periodicityConstraint,zeros(model.nStates+model.nControls,1),zeros(model.nStates+model.nControls,1),isSymmetric)
problem.addConstraint(@dynamicConstraintsRK,repmat(model.constraints.fmin,size(A,1)+1,N-1),repmat(model.constraints.fmax,size(A,1)+1,N-1), discMeth)
problem.addConstraint(@dynamicsFirstNodeConstraint,model.constraints.fmin,model.constraints.fmax)
end
