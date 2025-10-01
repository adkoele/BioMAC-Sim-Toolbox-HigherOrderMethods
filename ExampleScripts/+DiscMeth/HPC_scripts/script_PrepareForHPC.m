%======================================================================
%> @file +DiscMeth/script_PrepareForHPC.m
%> @brief Script to read osim files and adapt the mat files for HPC
%>
%> @details
%> This version overrights the old .mat files. It is therefore not
%> necessary to copy them individually to the HPC.
%>
%> @author Alexander Weiss adapted from Marlies Nitschke
%> @date July, 2025
%======================================================================

clear all; close all; clc;

%% Settings
% Where is the data?
workDirectoryHPC = '/home/hpc/iwb8/iwb8103h/DiscMeth_Workdirectory';
dataFolder      = '/data/DiscMeth';

% What to do?
createModel = 1;
adaptPath = 1;

%% Initialization
% Get path of this script
pathOfScript = mfilename('fullpath');
pathOfScript = strrep(pathOfScript, '\', '/');
pathOfScript = strsplit(pathOfScript, '/'); 

% Get path of the repository
path2repo = pathOfScript(1:end-3);         % use all but not src/scripts/+IMUTracking3D/script_PrepareModelForHPC
path2repo = strjoin(path2repo, filesep);   % connect cells

% Add path to repo to folders specified in the settings
dataFolderHPC      = [workDirectoryHPC dataFolder];
dataFolderPC       = [path2repo dataFolder];


%% Go over all model files
% Get all files
filelist = dir([dataFolderPC filesep '*' filesep '*.osim']); 

% Go over the files
for iFile = 1 : length(filelist)
    
    modelFile = [filelist(iFile).folder filesep filelist(iFile).name];
    
    if createModel
        % Create model (read osim and get moment arms)
        model = Gait3d(modelFile);
    end
    
    if adaptPath
        % Load old .mat file
        matFile = strrep(modelFile, '.osim', '.mat');
        load(matFile);
        
        % Adapt file name
        modelFileHPC = strrep(modelFile, dataFolderPC, dataFolderHPC);
        model.osim.file = modelFileHPC;
        
        % Save it
        save(matFile,'model','osim_sha256');
    end
    
end
