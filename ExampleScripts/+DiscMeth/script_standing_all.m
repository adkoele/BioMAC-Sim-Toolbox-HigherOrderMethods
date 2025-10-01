% ======================================================================
%> @file    DiscMeth/script_standing_all.m
%> @brief   Function to run standing simulations for a given participant.
%>
%> @details
%> This function runs the standing simulation for a specific participant 
%> using the provided trial data. The simulation can be performed in either 2D 
%> or 3D, depending on the 'dim' parameter.  
%>
%> @param   path2repo    String: Path to the repository 
%> @param   participant  String: Identifier for the participant whose data will be simulated
%> @param   trial_Npose  String: Trial name for the participant's static pose
%> @param   dim          String: Simulation dimension, either '2D' or '3D'
% ======================================================================

function script_standing_all(path2repo, participant, trial_Npose, dim)

    %% Run standing simulations
    % Standing with marker tracking
    resultFileStandingMarker = DiscMeth.run_standing_marker(path2repo, participant, trial_Npose, dim);
    
    % Standing for joint angle tracking
    resultFileStandingAngle = DiscMeth.run_standing_angle(path2repo, participant, trial_Npose, dim);

end