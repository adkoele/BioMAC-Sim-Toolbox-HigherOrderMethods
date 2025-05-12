


function script_standing_all(path2repo, participant, trial_Npose)

    % Error measure for analysis of error between tracking and simulated data
    measureNames = {'RMSE'};

    %% Run standing simulations
    % Standing with marker tracking
    resultFileStandingMarker = DiscMeth.run_standing_marker(path2repo, participant, trial_Npose);
    
    % Standing for joint angle tracking
    resultFileStandingAngle = DiscMeth.run_standing_angle(path2repo, participant, trial_Npose);
    
    % Compare standing simulations (can be easily rerun to inspect the results)
    % Load the marker tracking result
    load(resultFileStandingMarker);
    
    % Plot the stick figure visualizing the result
    figure();
    x = result.X(result.problem.idx.states);
    result.problem.model.showStick(x);
    markerTable = result.problem.objectiveTerms(1).varargin{1}.variables;
    markerMean = cell2mat(markerTable.mean')/1000; % extract and convert from mm to meter
    result.problem.model.showMarker(x, markerTable, markerMean);
    title('Standing Marker Tracking');
    view(45, 20);
    legend(findobj(gca, 'Type', 'Scatter'), {'Measured', 'Simulated'});
    
    % Extract simulated data
    settings.marker = markerTable;
    simVarTableStandingMarker = result.problem.extractData(result.X, settings);
    
    % Compute errors to tracked data
    simVarTableStandingMarker = calculateMeasures(simVarTableStandingMarker, measureNames);
    
    % Load the angle tracking result
    load(resultFileStandingAngle);
    
    % Plot the stick figure visualizing the result
    figure();
    x = result.X(result.problem.idx.states);
    result.problem.model.showStick(x);
    result.problem.model.showMarker(x, markerTable, markerMean);
    title('Standing Angle Tracking');
    view(45, 20);
    legend(findobj(gca, 'Type', 'Scatter'), {'Measured', 'Simulated'});
    
    % Extract simulated data
    simVarTableStandingAngle = result.problem.extractData(result.X, settings);
    
    % Compute errors to tracked data
    simVarTableStandingAngle = calculateMeasures(simVarTableStandingAngle, measureNames);

end