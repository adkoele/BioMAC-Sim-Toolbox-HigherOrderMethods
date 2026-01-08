% ======================================================================
%> @file    +DiscMeth/getError.m
%> @brief   Script to analyze discretization methods in 2D and 3D simulations.
%>
%> @details
%> This script evaluates the accuracy and performance of different 
%> discretization methods (e.g., BE, ME, LIIIc-2, RIIa-2, LIIIc-3, RIIa-3) 
%> across multiple participants, movements, and tracking strategies 
%> (Angle, Marker, IMU). It computes and compares dynamic residuals, 
%> joint errors, and timing information. The results are summarized 
%> in a table and can be visualized through error plots and summary plots.
%>
%>
%> @example
%>   % Run analysis in 2D for selected participants
%>   
%>
%>
%> @retval  dynResTable      Table summarizing dynamic residuals, errors, 
%>                           timing information, and convergence status 
%>                           for all tested discretization methods.
%>
%> @author  Alexander Weiss
%> @date    October, 2025
% ======================================================================

close all
clear all
clc

set(0, 'DefaultFigureWindowStyle', 'normal'); 



path2repo = what('BioMAC-Sim-Toolbox-HigherOrderMethods').path;
dataFolder = ['data', filesep, 'DiscMeth'];

resultFolder3D = ['results', filesep, 'DiscMeth3D'];

testDiscMethods = {'BE', 'ME','LIIIc-2', 'RIIa-2', 'LIIIc-3', 'RIIa-3'};

nPar = 11;
allPar = [2:nPar];
dim = '2D';

FontName = 'Times New Roman';
FontSize = 12*2;
LineWidth = 4;

% Interpolation points between collocation nodes
testJoint = 'knee_angle_r';%'pelvis_ty';%
plotError = 0;
plotAllTime = 0; % 0: plot Fraction, 1: plot First Order, 2: plot Second order, 3: plot third Order

plotAll = 1;

startNode = 0;
endNode = 10;

numInterp = 100;

dynResTable = table();
%column_width_cm = 10.5; % for single-column
column_width_cm = 19; % for two-column

% Convert cm to inches (MATLAB uses inches)
column_width_in = column_width_cm; %/ 2.54;
height_in = column_width_in;% * 0.75; % Aspect ratio (adjustable)

if strcmp(dim, '3D')
    resultFolder = ['results', filesep, 'DiscMeth3D'];
    movements = {'straightslowrunning', 'curvedslowrunning', 'vcut'};
    tracking = {'Angle', 'Marker', 'IMU'};
    nNodesList = [10, 25, 50, 75, 100];
    
elseif strcmp(dim, '2D')
    resultFolder = ['results', filesep, 'DiscMeth2D'];
    movements = {'straightslowrunning'};
    tracking = {'Angle', 'Marker', 'IMU'};
    nNodesList = [10, 25, 50, 75, 100, 150, 200];
end

if plotError == 1
    nPar = allPar(1);
    nNodesList = nNodesList(1);
    plotAll = 0;
    if plotAllTime == 1
        testDiscMethods = {'BE', 'ME'};
    elseif plotAllTime == 2
        testDiscMethods = {'LIIIc-2', 'RIIa-2'};
    elseif plotAllTime == 3
        testDiscMethods = {'LIIIc-3', 'RIIa-3'};
    end
end

%% Get Standard deviations
for iMov = 1%:size(movements,2)
    for iTra = 1%1:size(tracking,2)
        for iPar = 2:11
            if strcmp(dim, '3D')
                modelFile = [path2repo, filesep dataFolder,filesep,'Subject_',sprintf('%02d', iPar), filesep, 'Subject_', sprintf('%02d', iPar),'.osim'];
                model = Gait3d(modelFile);
            else
                modelFile2D = 'gait2dc_par.xls';
                model = Gait2dc(modelFile2D);
            end

            for iNodesTest = nNodesList(1)
                curResFolder = [path2repo,filesep,resultFolder,filesep,'Subject_',sprintf('%02d', iPar), filesep, movements{1,iMov},filesep, tracking{1,iTra}];
                for iDM = testDiscMethods
                    indDM = find(strcmp(iDM{1,1}, testDiscMethods));
                    try
                        resFileList(indDM,:) = dir(fullfile([curResFolder, filesep, iDM{1,1}, filesep, '*_', int2str(iNodesTest),'_nodes.mat']));
                    catch
                        continue
                    end
                end

                for iFile = 1:size(resFileList)
                    indDM = iFile;
                    try
                        resFile(iFile) = load([curResFolder, filesep, testDiscMethods{1, iFile}, filesep,resFileList(iFile,1).name]);
                    catch
                        continue
                    end
                end

                for iRes = 1:size(resFile,2)
                    resFileCur = resFile(iRes);%load(resFile(1,iRes).result.filename);
                    resultRunning = resFileCur.result;
                     X = resultRunning.X;
                    nNodes = resultRunning.problem.nNodesDur;
                    h = X(resultRunning.problem.idx.dur)/(nNodes-1);

                    settings.plotInitialGuess = 1;
                    style.figureSize = [0 0 16 26];
                    getFullCycle = 0;
                    dofNames = model.dofs.Properties.RowNames;
                    if strcmp(dim, '3D')
                        idxTrans = ismember(dofNames,{'pelvis_tx','pelvis_ty','pelvis_tz'});
                    else
                        idxTrans = ismember(dofNames,{'pelvis_tx','pelvis_ty'});
                    end
                    dofAngles = dofNames(~idxTrans);
                    settings.angle = dofAngles;
                    simVarTable = resultRunning.problem.extractData(resultRunning.X, settings, [], getFullCycle);
                    meanTestJoint(iPar-1,iRes,:,:) = simVarTable.mean(strcmp(simVarTable.type, 'angle')); %(strcmp(simVarTable.name, testJoint) & 
                    TESTmeanTestJoint(iPar-1,iRes,:,:) = simVarTable.sim(strcmp(simVarTable.type, 'angle')); %(strcmp(simVarTable.name, testJoint) & 

                    for k = 1:size(meanTestJoint,3)
                        % Participants, Methods, Joints, Nodes
                        avgData(iPar-1,iRes,k,:) = meanTestJoint{iPar-1,iRes,k};
                        TESTavgData(iPar-1,iRes,k,:) = TESTmeanTestJoint{iPar-1,iRes,k};

                    end
                end
            end
        end
    end
end

% Derivative of tracking data
avgDatadt = diff(avgData,1,4)./h;
avgDatadt = avgDatadt./180.*pi;
avgDatadtdt = diff(avgDatadt,1,4)./h;

% Mean and std of tracking data over participants
mean_avgDatadt = mean(avgDatadt,1);
std_avgDatadt = std(avgDatadt,0,1);
meanVals = reshape(mean_avgDatadt, size(mean_avgDatadt,2), size(mean_avgDatadt,3), size(mean_avgDatadt,4));
stdVals = reshape(std_avgDatadt, size(std_avgDatadt,2), size(std_avgDatadt,3), size(std_avgDatadt,4));
% Mean of stds over methods (should all be equal)
stdSumVals = squeeze(mean(stdVals,1));

meanValsdt = squeeze(mean(avgDatadtdt,1));
stdValsdt = squeeze(std(avgDatadtdt,0,1));
stdSumValsdt = squeeze(mean(stdVals,1));

% Sum of stds over joints 
sumJointsstd = sum(abs(stdSumVals),1);% / size(stdSumVals,1);
intSumJoints = sum(sumJointsstd) / size(sumJointsstd,2);
%end

for iMov = 1:size(movements,2)
    for iTra = 1%:size(tracking,2)
        for iPar = 2:nPar
            clear iNodesTest

            if strcmp(dim, '3D')
                modelFile = [path2repo, filesep dataFolder,filesep,'Subject_',sprintf('%02d', iPar), filesep, 'Subject_', sprintf('%02d', iPar),'.osim'];
                model = Gait3d(modelFile);
            else
                modelFile2D = 'gait2dc_par.xls';
                model = Gait2dc(modelFile2D);
            end

            for iNodesTest = nNodesList
                clear x1_int_T x2_int_T xc_int_T u_int_T u2_int_T
                clear errorTimeTrans errorTimeRad errorMethRad errorMethTrans errorSumRad errorSumTrans errorSumRaddot errorSumTransdot...
                    intErrorRad intErrorTrans intErrorRaddot intErrorTrans meanTestJoint
                curResFolder = [path2repo,filesep,resultFolder,filesep,'Subject_',sprintf('%02d', iPar), filesep, movements{1,iMov},filesep, tracking{1,iTra}];
                for iDM = testDiscMethods
                    indDM = find(strcmp(iDM{1,1}, testDiscMethods)); 
                    try
                        resFileList(indDM,:) = dir(fullfile([curResFolder, filesep, iDM{1,1}, filesep, '*_', int2str(iNodesTest),'_nodes.mat']));
                    catch
                        continue
                    end
                end

                for iFile = 1:size(resFileList)
                    indDM = iFile;
                    try
                        resFile(iFile) = load([curResFolder, filesep, testDiscMethods{1, iFile}, filesep,resFileList(iFile,1).name]);
                    catch
                        continue
                    end
                end

                for iRes = 1:size(resFile,2)
                    resFileCur = resFile(iRes);%load(resFile(1,iRes).result.filename);
                    resultRunning = resFileCur.result;
                    settings.plotInitialGuess = 1;
                    style.figureSize = [0 0 16 26];
                    getFullCycle = 0;
                    dofNames = model.dofs.Properties.RowNames;
                    if strcmp(dim, '3D')
                        idxTrans = ismember(dofNames,{'pelvis_tx','pelvis_ty','pelvis_tz'});
                    else
                        idxTrans = ismember(dofNames,{'pelvis_tx','pelvis_ty'});
                    end

                    dofAngles = dofNames(~idxTrans);
                    settings.angle = dofAngles;
                    simVarTable = resultRunning.problem.extractData(resultRunning.X, settings, [], getFullCycle);
                    indJoint = find(strcmp(resultRunning.problem.model.states.name, testJoint) & strcmp(resultRunning.problem.model.states.type, 'q'));
                    simVarTableJoints = simVarTable(strcmp(simVarTable.type, 'angle'),:);
                    indJointSim = find(strcmp(simVarTableJoints.name, testJoint));%  & (strcmp(simVarTable.type, 'translation') | strcmp(simVarTable.type, 'angle')));
                    iCurNodes = find(nNodesList==iNodesTest);
                   
                        
                    X = resultRunning.X;
                    nNodes = resultRunning.problem.nNodesDur;
                    h = X(resultRunning.problem.idx.dur)/(nNodes-1);
     
                    idx_q = 1:resultRunning.problem.model.nDofs;
                    idx_qdot = resultRunning.problem.model.nDofs+1:2*resultRunning.problem.model.nDofs;
                    cp_ind = find(strcmp(resultRunning.problem.model.constraints.type, 'CP'));


                    dynResAll = array2table(zeros(size(testDiscMethods,2),21));
                    dynResAll.Properties.VariableNames = {'Method','DynRes', 'DynResDot', 'nNodes', ...
                        'ErrorT', 'ErrorR', 'SErrorR', 'SErrorT', 'WallTime', 'CpuTime',...
                        'Participant', 'Movement', 'Tracking',...
                        'SumErrorR', 'SumErrorT', 'SumSecErrorR', 'SumSecErrorT',...
                        'IntErrorR', 'IntErrorT', 'IntSecErrorR', 'IntSecErrorT'};
                    dynResAll.Properties.RowNames = {};

                    [x_all_interp, xdot_all_interp, residual_all_interp, u_all_interp] = DiscMeth.interpolateRK(resultRunning, numInterp);

                    indStat = strcmp(resFile(iRes).result.problem.model.states.type, 'q');
                    indStatdot = strcmp(resFile(iRes).result.problem.model.states.type, 'qdot');
                    indS = strcmp(resFile(iRes).result.problem.model.states.type, 's');
                    indA = strcmp(resFile(iRes).result.problem.model.states.type, 'a');

                    indTrans = contains(resFile(iRes).result.problem.model.states.name, 'pelvis_t') & ...
                        ~contains(resFile(iRes).result.problem.model.states.name, 'pelvis_tilt') & ...
                        strcmp(resFile(iRes).result.problem.model.states.type, 'q');

                    discMethAll{iRes} = resultRunning.problem.Euler;
                    dynResVal{iRes} = residual_all_interp(idx_q,:)';
                    dynResValdot{iRes} = residual_all_interp(idx_qdot,:)';
                    nNodesAll(iRes,:) = nNodes;
                    errorTimeTrans(iRes,:) = sum(abs(residual_all_interp(1,:)),2);
                    errorTimeRad(iRes,:) = sum(abs(residual_all_interp(indJoint,:)),2); %!!!!

                    errorTimeAct(iRes,:) = sum(abs(residual_all_interp(indA,:)),2);
                    errorTimeCEL(iRes,:) = sum(abs(residual_all_interp(indS,:)),2); %!!!!

                    errorMethTrans(iRes,:) = sum(errorTimeTrans(iRes,:)) / size(residual_all_interp,2);
                    errorMethRad(iRes,:) = sum(errorTimeRad(iRes,:)) / size(residual_all_interp,2);

                    errorMethAct(iRes,:) = sum(errorTimeAct(iRes,:)) / size(residual_all_interp,2);
                    errorMethCEL(iRes,:) = sum(errorTimeCEL(iRes,:)) / size(residual_all_interp,2);
                    wallTime(iRes) = resultRunning.info.wallTime;
                    cpuTime(iRes) = resultRunning.info.cpu;
                    convergedAll(iRes) = resultRunning.converged;

                    

                    resRad = residual_all_interp(idx_q,:);
                    resTrans = resRad(indTrans,:);
                    resRad(indTrans,:) = [];

                    resRaddot = residual_all_interp(idx_qdot,:);
                    resTransdot = resRaddot(indTrans,:);
                    resRaddot(indTrans,:) = [];

                    resS = residual_all_interp(indS,:);
                    resA = residual_all_interp(indA,:);
                    
                    errorSumRad{iRes} = sum(abs(resRad),1);
                    errorSumTrans{iRes} = sum(abs(resTrans),1);
                    errorSumS{iRes} = sum(abs(resS),1);
                    errorSumA{iRes} = sum(abs(resA),1);

                    intErrorRad(iRes) = sum(errorSumRad{iRes}) / size(residual_all_interp,2);
                    intErrorTrans(iRes) = sum(errorSumTrans{iRes}) / size(residual_all_interp,2);

                    intErrorA(iRes) = sum(errorSumA{iRes}) / size(residual_all_interp,2);
                    intErrorS(iRes) = sum(errorSumS{iRes}) / size(residual_all_interp,2);
                    
                    errorSumRaddot{iRes} = sum(abs(resRaddot),1);
                    errorSumTransdot{iRes} = sum(abs(resTransdot),1);

                    intErrorRaddot(iRes) = sum(errorSumRaddot{iRes}) / size(residual_all_interp,2);
                    intErrorTransdot(iRes) = sum(errorSumTransdot{iRes}) / size(residual_all_interp,2); 

                    hAll(iRes) = h;

                    secerrorTimeTrans(iRes,:) = sum(abs(residual_all_interp(1+size(idx_q,2),:)),2);
                    secerrorTimeRad(iRes,:) = sum(abs(residual_all_interp(indJoint+size(idx_q,2),:)),2); %!!!!
                    secerrorMethTrans(iRes,:) = sum(secerrorTimeTrans(iRes,:)) / size(residual_all_interp,2);
                    secerrorMethRad(iRes,:) = sum(secerrorTimeRad(iRes,:)) / size(residual_all_interp,2);

                    iMovement{iRes} = movements{1,iMov};
                    iTracking{iRes} = tracking{1,iTra};

                    meanTestJoint(iRes,:) = simVarTable.mean((strcmp(simVarTable.name, testJoint) & strcmp(simVarTable.type, 'angle')));
                    intThres(iRes) = intSumJoints;

                    %errorMethAll{iRes} = sum(abs(errorMeth(iRes,:))) / size(errorMeth,2);
                end

                % Specification
                dynResAll.Method = discMethAll';
                dynResAll.nNodes = nNodesAll;
                dynResAll.Movement = iMovement';
                dynResAll.Tracking = iTracking';
                dynResAll.Participant = iPar*ones(size(discMethAll,2),1);
                % Dynamic residuals
                dynResAll.DynRes = dynResVal';
                dynResAll.DynResDot = dynResValdot';
                % Errors one joint
                dynResAll.ErrorA = errorMethAct;
                dynResAll.ErrorS = errorMethCEL;
                dynResAll.ErrorT = errorMethTrans;
                dynResAll.ErrorR = errorMethRad;
                dynResAll.SErrorR = secerrorMethRad;
                dynResAll.SErrorT = secerrorMethTrans; 
                % Sum Error
                dynResAll.SumErrorR = errorSumRaddot';
                dynResAll.SumErrorT = errorSumTransdot';
                dynResAll.SumSecErrorR = errorSumRad';
                dynResAll.SumSecErrorT = errorSumTrans';
                dynResAll.IntErrorR = intErrorRad';
                dynResAll.IntErrorT = intErrorTrans';
                dynResAll.IntSecErrorR = intErrorRaddot';
                dynResAll.IntSecErrorT = intErrorTransdot';
                % Muscle Errors
                dynResAll.SumErrorA = errorSumA';
                dynResAll.SumErrorS = errorSumS';
                dynResAll.IntErrorA = intErrorA';
                dynResAll.IntErrorS = intErrorS';
                % Timing
                dynResAll.WallTime = wallTime';
                dynResAll.CpuTime = cpuTime';
                dynResAll.Converged = convergedAll';
                % Tracking Data
                dynResAll.tracking = meanTestJoint;
                dynResAll.deltaT = hAll';
                dynResAll.IntThreshold = intThres';

                dynResTable = [dynResTable; dynResAll];

                for iMeth = 1:size(dynResAll,1)
                    disp([dynResAll.Method{iMeth,1}, '_rad: ' num2str(dynResAll.ErrorR(iMeth)), ' for Nodes: ',  int2str(iNodesTest)]);
                    disp([dynResAll.Method{iMeth,1}, '_trans: ' num2str(dynResAll.ErrorT(iMeth)), ' for Nodes: ',  int2str(iNodesTest)]);
                end

                if plotError == 1
                    DiscMeth.plotError(curResFolder, resultRunning, dynResAll, column_width_in, height_in, numInterp, indJoint, startNode, endNode, FontName, FontSize, LineWidth, plotAllTime, stdSumVals,stdSumValsdt, indJointSim, movements{1,iMov})
                end
            end            
        end
        %DiscMeth.analyzeANOVA([path2repo, filesep, resultFolder], dynResTable, movements{1,iMov}, dim,FontName, FontSize);

        
        if plotAll == 1
            DiscMeth.plotAll([path2repo,filesep,resultFolder], dynResTable, column_width_in, height_in, nPar,movements{1,iMov}, tracking{1,iTra}, FontName, FontSize, LineWidth, dim)
        end
    end
end





