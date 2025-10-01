%======================================================================
%> @file +DiscMeth/run_standing_marker.m
%> @brief Function to do standing simulation with marker tracking
%>
%> @author Alexander Weiss adapted from Marlies Nitschke
%> @date December, 2021, July 2025
%======================================================================

%======================================================================
%> @brief Function to do standing simulation with marker tracking
%>
%> @details
%> This function tracks data during standing to optimize the CP position.
%>
%> @param  workDirectory    String: Current work directory
%> @param  participant      String: Name of participant (e.g. Participant_02)
%> @param  trial            String: Name of trial (e.g. trial0002)
%> @retval resultFile       String: Filename of result
%======================================================================
function resultFile = run_standing_marker(workDirectory, participant, trial, dim)

% Fixed settings
dataFolder     = ['data' filesep 'DiscMeth']; % Relative from the work directory
dataTrackFile  = [participant filesep 'N-Pose' filesep trial '_DataStructMeasured.mat'];
if strcmp(dim, '2D')
    modelFile   = 'gait2dc_par.xls'; 
    resultFolder   = ['results' filesep 'DiscMeth2D' filesep participant];  % Relative from the path of the repository
elseif strcmp(dim, '3D')
    modelFile      = [participant filesep participant '.osim'];
    resultFolder   = ['results' filesep 'DiscMeth3D' filesep participant];  % Relative from the path of the repository
else
    error('Please specify either 2D or 3D as dimension');
end

% Get absolute file names
resultFile     = [workDirectory filesep resultFolder filesep 'standing_marker'];
dataTrackFile  = [workDirectory filesep dataFolder filesep dataTrackFile];
%modelFile      = [workDirectory filesep dataFolder filesep modelFile];
if strcmp(dim, '3D')
    modelFile    = [workDirectory filesep dataFolder filesep modelFile];
end
% Create resultfolder if it does not exist
if ~exist([workDirectory,filesep,resultFolder], 'dir')
    mkdir([workDirectory,filesep,resultFolder]);
end

% Create an instane of our 2D or 3D model class
if strcmp(dim, '2D')
    model = Gait2dc(modelFile);
elseif strcmp(dim, '3D')
    model = Gait3d(modelFile);
else
    error('Please specify either 2D or 3D as dimension');
end

% Solve simulation for multiple random initial guesses
rng('default');
nRep = 10;
results = cell(nRep, 1);
objSum = nan(nRep, 1);
for iRep = 1 : nRep
    % Specify the optimizaton problem
    iSample = 10; % Index of the used sample in the data
    W.trackMarker= 1e-2;        % Weight of marker tracking term in objective
    W.trackGRF   = 1e-2;        % Weight of GRF tracking term in objective
    W.effMuscles = 1e+00;       % Weight of effort term for muscles in objective
    W.effTorques = 1e-01;       % Weight of effort term for torques in objective
    problemStanding = DiscMeth.setup_standing_marker(model, dataTrackFile, iSample, resultFile, W);
   % problemStanding = DiscMeth.standing_2D(model, resultFile);
    % Solve
    solver = IPOPT();
    solver.setOptionField('max_iter', 20000);
    solver.setOptionField('tol', 0.0001);
    results{iRep} = solver.solve(problemStanding);
    objSum(iRep) = sum([results{iRep}.problem.objectiveTerms.weightedValue]);
end

% Save simulation result of result with lowest objective value
[~, iMinObj] = min(objSum);
result = results{iMinObj};
result.save(resultFile);
end
