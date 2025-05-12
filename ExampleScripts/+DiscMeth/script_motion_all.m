%======================================================================
%> @file MarkerTracking3D/scriptMarkerTracking.m
%> @brief Example script for marker tracking (in 3D)
%> @example
%> @details
%> This is a script that shows how to generate marker tracking simulations
%> in 3D without any task constraints. These scripts are based on the
%> following paper: Nitschke, M., Marzilger, R., Leyendecker, S., Eskofier,
%> B. M., & Koelewijn, A. D. (2023). Change the direction: 3D optimal 
%> control simulation by directly tracking marker and ground reaction 
%> force data. PeerJ, 11, e14852.
%> Additional data can be found here: DOI 10.5281/zenodo.6949011. We only
%> use one trial of each condition for participant 2.
%>
%> @author Marlies Nitschke, Anne Koelewijn
%> @date September, 2024
%======================================================================

clear all 
close all
clc

%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));

% Path to your repository
path2repo = [filePath filesep '..' filesep '..'];


% Loop through all data

%% Settings
%workDirectory = '/home/qe21vady/Documents/Code/projects-simulation-project3d';
% Range can be limited to specific subjects
rangePar = 3:11;%2:11; % Subject 1 has missing marker
% Range can be limited to specific movements
rangeMov = [2];%, 4, 5];
% Range can be limited to specific Discretization methods
rangeDisc = 1:2;%6;
% Range can be limted to specific number of Nodes
rangeNodes = [15, 25];%, 50, 75, 100];%, 150, 200, 250, 400];

% Weights
WMar = 10^-2; %logspace(-3, -1, 3);
WGRF = 10^-3; %logspace(-3, -1, 3);
WReg = 1e-3;

%% Initialization
% Data name of the ten subjects
participantCell = DiscMeth.getConfigStudyParameters('participantCell');
% Trial name of NPose
staticCell = DiscMeth.getConfigStudyParameters('NPoseTrialNamesCell');
% Movement types used in names of files
movementTypeCell = DiscMeth.getConfigStudyParameters('movementTypesCell');
% Trial names for each movement and subject
trialNamesCell = DiscMeth.getConfigStudyParameters('OMCTrialNamesCell');
% Discretization methods names
discMethCell = DiscMeth.getConfigStudyParameters('discMethCell');
% Make weights to strings for input (on HPC we need strings)
WMar = cellfun(@(x) sprintf('%.0e', x), num2cell(WMar), 'UniformOutput', false);
WGRF = cellfun(@(x) sprintf('%.0e', x), num2cell(WGRF), 'UniformOutput', false);
WReg = cellfun(@(x) sprintf('%.0e', x), num2cell(WReg), 'UniformOutput', false);

% Function to call Butcher Tableau
% Define number of collocation Nodes here



% participant = 'Participant_02';
% % Trial names
% trial_Npose = 'trial0002';
% trial_straightrunning = 'trial0025';
% trial_curvedslowrunning = 'trial0100';
% trial_vcut = 'trial0122';
% movements = {'vcut'};%'straightrunning', 'curvedrunning', 'vcut'};
% trials = {trial_vcut};%trial_straightrunning, trial_curvedslowrunning, trial_vcut};



% Error measure for analysis of error between tracking and simulated data
measureNames = {'RMSE'};

% Dock all figures in one Window
set(0, 'DefaultFigureWindowStyle', 'docked'); 

% Adapt style (You might have to adapt this to best fit your screen)
style.subFigSettings.nCol = 12;
style.subFigSettings.width = 3;
style.subFigSettings.height = 2;
style.trackColor = 'k--';
style.extraColor = 'r:';
style.xLabelText  = 'Motion in \%';
style.lineWidth = 1.5;


for iPar = rangePar

    participant = participantCell{iPar};
    trial_Npose = staticCell{iPar};
    %% Run standing simulation
    DiscMeth.script_standing_all(path2repo, participant, trial_Npose)

    %% Run movement simulations
    % Simulate the different movements as specified before. This will require a
    % couple of hours per simulation.
    for iMov = rangeMov
        movement = movementTypeCell{iMov};
        trials = trialNamesCell{iPar, iMov};
        for iDisc = rangeDisc
            discMeth = discMethCell{iDisc};
            for iNod = 1:numel(rangeNodes)
                nNodes = rangeNodes(iNod);
                % Marker tracking
                resultFilesRunningMarker{iPar, iMov, iNod, iDisc} = DiscMeth.run_motion_marker(path2repo, participant, movement, trials(1,:), discMeth, nNodes);
                % Joint angle tracking
                resultFilesRunningAngle{iPar, iMov, iNod, iDisc} = DiscMeth.run_motion_angle(path2repo, participant, movement, trials(1,:), discMeth, nNodes);
            end
        end
    end
end

%% Compare movement simulations (can be easily rerun to inspect the results)
% Get simulated variables
% for iMov = 1:length(movements)
% 
%     % Load the marker tracking result
%     load(resultFilesRunningMarker{iMov});
% 
%     % Extract simulated data
%     markerTable = result.problem.objectiveTerms(1).varargin{1}.variables;
%     settings.marker = markerTable;
%     simVarTableRunningMarker{iMov} = result.problem.extractData(result.X, settings);
% 
%     % Compute errors to tracked data
%     simVarTableRunningMarker{iMov} = calculateMeasures(simVarTableRunningMarker{iMov}, measureNames);
% 
% 
%     % Load the angle tracking result
%     load(resultFilesRunningAngle{iMov});
% 
%     % Extract simulated data
%     simVarTableRunningAngle{iMov} = result.problem.extractData(result.X, settings);
% 
%     % Compute errors to tracked data
%     simVarTableRunningAngle{iMov} = calculateMeasures(simVarTableRunningAngle{iMov}, measureNames);
% 
%     % Make plots
%     Collocation.plotMultSimVarTables({simVarTableRunningMarker{iMov}, simVarTableRunningAngle{iMov}}, style);
%     legend({'Marker', 'Angle'})
% 
%     % Pause loop to look at results
%     if iMov < length(movements)
%         fprintf('Inspect the results for %s. Press a key to close the figures and continue with the next movement. \n', movements{iMov})
%         pause;
%         close all;
%     end
% end



