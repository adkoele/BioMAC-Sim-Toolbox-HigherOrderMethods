%======================================================================
%> @file    DiscMeth/script_motion_all.m
%> @brief   Script to run 2D or 3D simulations over multiple movements 
%>          and tracking strategies to test discretization methods.
%>
%> @details
%> The script loops over a configurable set of participants, movement 
%> types, discretization methods, and numbers of nodes. It runs both 
%> standing simulations and dynamic movement simulations (marker tracking 
%> and joint angle tracking) in 2D or 3D.  
%>  
%> The primary purpose is to evaluate discretization schemes (e.g. BE, ME, 
%> LIIIc-2, RIIa-2, LIIIC-3, RIIa-3) under different conditions.  
%>  
%> @author  Alexander Weiss
%> @date    July 2025
%======================================================================

clear all 
close all
clc

%% Settings
% Get path of this script
filePath = fileparts(mfilename('fullpath'));

% Path to your repository
path2repo = what('BioMAC-Sim-Toolbox-HigherOrderMethods').path;

%Define dimension to be simulated
dim = '3D';

% Loop through all data
%% Settings
% Range can be limited to specific subjects
rangePar = [2:11];
% Range can be limited to specific movements
rangeMov = [2, 4, 5];
% Range can be limited to specific Discretization methods
rangeDisc = [1:6]; % 1:BE, 2:ME, 3:LIIIc-2, 4:RIIa-2, 5:LIIIC-3, RIIa-3
% Range can be limted to specific number of Nodes
rangeNodes = [10, 25, 50, 75, 100]; %, 150, 200];

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

% Iterate over participants. Potentially parfor
for idxPar = 1:numel(rangePar)
   %try
    iPar = rangePar(idxPar);
    participant = participantCell{iPar};
    trial_Npose = staticCell{iPar};
    %% Run standing simulation
    DiscMeth.script_standing_all(path2repo, participant, trial_Npose, dim)

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
                % Marker tracking -> Not used in the paper
                 if strcmp(dim, '2D')
                    DiscMeth.run_motion_marker_2D(path2repo, participant, movement, trials(1,:), discMeth, nNodes);
                 elseif strcmp(dim, '3D')
                    DiscMeth.run_motion_marker_3D(path2repo, participant, movement, trials(1,:), discMeth, nNodes);
                 end
                % Joint angle tracking 
                if strcmp(dim, '2D')
                    DiscMeth.run_motion_angle_2D(path2repo, participant, movement, trials(1,:), discMeth, nNodes);
                elseif strcmp(dim, '3D')
                    DiscMeth.run_motion_angle_3D(path2repo, participant, movement, trials(1,:), discMeth, nNodes);
                end
            end
        end
    end
end






