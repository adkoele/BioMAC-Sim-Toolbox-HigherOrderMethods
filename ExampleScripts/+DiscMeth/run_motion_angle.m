%======================================================================
%> @file +MarkerTracking3D/run_motion_angle.m
%> @brief Function to do simulation with angle tracking
%>
%> @author Marlies Nitschke
%> @date December, 2021
%======================================================================

%======================================================================
%> @brief Function to do standing simulation with angle tracking
%>
%> @details
%> This function tracks data to simulate a motion without periodicity.
%>
%> @param  workDirectory    String: Current work directory
%> @param  participant          String: Name of participant (e.g. Participant_02)
%> @param  motion           String: Name of motion (e.g. straightslowrunning)
%> @param  trial            String: Name of triak (e.g. trial0025)
%> @retval resultFile       String: Filename of result
%======================================================================
function resultFile = run_motion_angle(workDirectory, participant, motion, trial, discMeth, nNodes)

% Fixed settings
dataFolder     = ['data' filesep 'DiscMeth']; % Relative from the work directory
dataMeasFile   = [participant filesep motion filesep trial '_DataStructMeasured.mat'];
dataInvFile    = [participant filesep motion filesep trial '_DataStructInverse.mat'];
modelFile      = [participant filesep participant '.osim'];
resultFolderStanding   = ['results' filesep 'DiscMeth3D' filesep participant];  % Relative from the path of the repository
resultFolder   = ['results' filesep 'DiscMeth3D' filesep participant filesep motion filesep 'Angle' filesep discMeth];  % Relative from the path of the repository
resultFileStanding = 'standing_angle';
resultFile     = sprintf('%s_angle_%s_%s_nodes', trial, discMeth, int2str(nNodes));

% Get absolute file names
resultFileStanding = [workDirectory filesep resultFolderStanding filesep resultFileStanding];
resultFile         = [workDirectory filesep resultFolder filesep resultFile];
dataMeasFile       = [workDirectory filesep dataFolder filesep dataMeasFile];
dataInvFile        = [workDirectory filesep dataFolder filesep dataInvFile];
modelFile          = [workDirectory filesep dataFolder filesep modelFile];

% Create resultfolder if it does not exist
if ~exist([workDirectory,filesep,resultFolder], 'dir')
    mkdir([workDirectory,filesep,resultFolder]);
end

% Create an instane of our 3D model class
model = Gait3d(modelFile);

% Adjust CPs based on standing solution
resultStanding = load(resultFileStanding);
model.CPs.position = model.CPs.position + repmat([0 resultStanding.result.X(resultStanding.result.problem.idx.CPYOffset) 0], model.nCPs, 1);

% Specify the optimizaton problem
W.trackTrans = 1e-3;         % Weight of translation tracking term in objective
W.trackAngle = 1e-1;         % Weight of angle tracking term in objective
W.trackGRF   = 1e-3;         % Weight of GRF tracking term in objective
W.effMuscles = 1e+00;        % Weight of effort term for muscles in objective
W.effTorques = 1e-01;        % Weight of effort term for torques in objective
W.reg        = 1e-3;         % Weight of regularization term in objective
initialGuess = resultFileStanding;
nNodesEarlier = 10; %Number of time points before the actual motion of interest to ensure motion artefacts are before the motion of interest
problemRunning = DiscMeth.setup_motion_angle(model,dataInvFile,dataMeasFile,initialGuess,resultFile,W, nNodesEarlier, discMeth, nNodes);

% Solve
solver = IPOPT();
solver.setOptionField('max_iter', 2);
solver.setOptionField('tol', 0.0001);
resultRunning = solver.solve(problemRunning);
resultRunning.save(resultFile);

end
