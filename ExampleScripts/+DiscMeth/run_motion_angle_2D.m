%======================================================================
%> @file +DiscMeth/run_motion_angle_2D.m
%> @brief Function to do simulation with angle tracking
%>
%> @author Alexander Weiss adapted from Marlies Nitschke
%> @date July, 2025
%======================================================================

%======================================================================
%> @brief Function to do standing simulation with angle tracking
%>
%> @details
%> This function tracks data to simulate a motion without periodicity.
%>
%> @param  workDirectory    String: Current work directory
%> @param  participant      String: Name of participant (e.g. Participant_02)
%> @param  motion           String: Name of motion (e.g. straightslowrunning)
%> @param  trial            String: Name of triak (e.g. trial0025)
%> @param  discMeth         String: Name of discretization method (e.g. BE or RIIa-3)
%> @param  nNodes           Int or String: Number of Nodes (e.g. 50 or '50')
%> @retval resultFile       String: Filename of result
%======================================================================
function resultFile = run_motion_angle_2D(workDirectory, participant, motion, trial, discMeth, nNodes)

% Fixed settings
dataFolder     = ['data' filesep 'DiscMeth']; % Relative from the work directory
dataMeasFile   = [participant filesep motion filesep trial '_DataStructMeasured.mat'];
dataInvFile    = [participant filesep motion filesep trial '_DataStructInverse.mat'];
modelFile2D    = 'gait2dc_par.xls';
resultFolderStanding   = ['results' filesep 'DiscMeth2D' filesep participant];  % Relative from the path of the repository
resultFolder   = ['results' filesep 'DiscMeth2D' filesep participant filesep motion filesep 'Angle' filesep discMeth];  % Relative from the path of the repository
resultFileStanding = 'standing_angle_2D';
if isa(nNodes, 'char')
    resultFile     = sprintf('%s_angle_%s_%s_nodes', trial, discMeth, nNodes);
else
    resultFile     = sprintf('%s_angle_%s_%s_nodes', trial, discMeth, int2str(nNodes));
end

% Get absolute file names
resultFileStanding = [workDirectory filesep resultFolderStanding filesep resultFileStanding];
resultFile         = [workDirectory filesep resultFolder filesep resultFile];
dataMeasFile       = [workDirectory filesep dataFolder filesep dataMeasFile];
dataInvFile        = [workDirectory filesep dataFolder filesep dataInvFile];
%modelFile          = [workDirectory filesep dataFolder filesep modelFile];

% Create resultfolder if it does not exist
if ~exist([workDirectory,filesep,resultFolder], 'dir')
    mkdir([workDirectory,filesep,resultFolder]);
end

% Create an instane of our 2D model class
model = Gait2dc(modelFile2D);

% Adjust CPs based on standing solution
resultStanding = load(resultFileStanding);

% Specify the optimizaton problem
W.trackTrans = 1;         % Weight of translation tracking term in objective
W.trackAngle = 0.6;         % Weight of angle tracking term in objective
W.trackGRF   = 0.4;         % Weight of GRF tracking term in objective
W.effMuscles = 1e+02;        % Weight of effort term for muscles in objective
W.reg        = 1e-4;         % Weight of regularization term in objective
initialGuess = resultFileStanding;
nNodesEarlier = 10; %Number of time points before the actual motion of interest to ensure motion artefacts are before the motion of interest
if isa(nNodes,'char')
    nNodes = str2num(nNodes);
end

%Define angle tracking problem
problemRunning = DiscMeth.setup_motion_angle_2D(model,dataInvFile,dataMeasFile,initialGuess,resultFile,W, nNodesEarlier, discMeth, nNodes);

%Solve
solver = IPOPT();
solver.setOptionField('hsllib', '/usr/local/lib/libcoinhsl.dylib');
solver.setOptionField('linear_solver', 'ma57');
solver.setOptionField('max_iter', 20000);
solver.setOptionField('tol', 0.0001);
resultRunning = solver.solve(problemRunning);
resultRunning.save(resultFile);

end