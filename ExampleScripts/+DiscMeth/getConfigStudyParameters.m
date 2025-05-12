%======================================================================
%> @file +MarkerTracking3D/getConfigStudyParameters.m
%> @brief Get all study parameters which have to be entered here, such that
%> they don't have to be entered multiple times in the code.
%>
%> @todo If there is time. Think about a better structure saved in files.
%> Or add some error checking.
%>
%> @author Marlies Nitschke, Alexander Weiss
%> @date May, 2025
%======================================================================
%======================================================================
%> @brief Get all study parameters which have to be entered here
%>
%> @details
%> New function parameters are added to this file. Afterwards they can be
%> read easily from everywhere in the project if this file is in the MATLAB
%> path. Make sure if you enter new values, that e.g. subject{i}
%> corresponds always to weight{i}, and movementTypes{j} corresponds to
%> triggerIndices{j}.
%>
%> @param key to access the value. Describtions can be found in the code.
%> Possible values are:
%>      - subjectCell               (length = Number of subjects)
%>      - NPoseTrialNamesCell       (length = Number of subjects)
%>      - movementTypesCell         (length = Number of movement types)
%>      - movementNamesCell         (length = Number of movement types)
%>      - OMCTrialNamesCell         (size() = [Numb. sujects, Numb. movements])
%>
%> @retval value contains the data depending on the entered key
%======================================================================
function value = getConfigStudyParameters(key)
if strcmp(key, 'participantCell') 
    value = getParticipantCell();   
% cells with length = number of movements
elseif strcmp(key, 'movementTypesCell')
    value = getMovementTypesCell();
elseif strcmp(key, 'movementNamesCell')
    value = getMovementNamesCell();
elseif strcmp(key, 'discMethCell')
    value = getDiscMethCell();
elseif strcmp(key, 'NPoseTrialNamesCell')
    value = getNPoseTrialNamesCell();
elseif strcmp(key, 'OMCTrialNamesCell')
    value = getOMCTrialNamesCell();
% invalid key
else
    error('getConfigStudyParameters(): Invalid key was used');
end
end % end function getConfigStudyParameters()
%======================================================================
%> @brief Get all subject file names
%>
%> @retval value is cell array with name strings of the subject files
%======================================================================
function value = getParticipantCell()
    subjects = cell(11, 1);
    subjects{1} = 'Subject_01';  
    subjects{2} = 'Subject_02';   
    subjects{3} = 'Subject_03';  
    subjects{4} = 'Subject_04';  
    subjects{5} = 'Subject_05';  
    subjects{6} = 'Subject_06';  
    subjects{7} = 'Subject_07'; 
    subjects{8} = 'Subject_08';  
    subjects{9} = 'Subject_09';  
    subjects{10} ='Subject_10';
    subjects{11} ='Subject_11';
    value = subjects;
end
%======================================================================
%> @brief Get all movement types
%>
%> @retval value is cell array with movement types (file and path names)
%======================================================================
function value = getDiscMethCell()
    discMeth = cell(6, 1);
    
    discMeth{1} = 'BE';
    discMeth{2} = 'ME';
    discMeth{3} = 'LIIIc-2';
    discMeth{4} = 'RIIa-2';
    discMeth{5} = 'LIIIc-3';
    discMeth{6} = 'RIIa-3';
    
    value = discMeth;
end
%===
%======================================================================
%> @brief Get all movement types
%>
%> @retval value is cell array with movement types (file and path names)
%======================================================================
function value = getMovementTypesCell()
    movement = cell(5, 1);
    
    movement{1} = 'straightwalking';
    movement{2} = 'straightslowrunning';
    movement{3} = 'straightfastrunning';
    movement{4} = 'curvedslowrunning';
    movement{5} = 'vcut';
    
    value = movement;
end
%======================================================================
%> @brief Get all movement names
%>
%> @retval value is cell array with movement names (nice names for plotting)
%======================================================================
function value = getMovementNamesCell()
    movement = cell(5, 1);
    
    movement{1} = 'Straight Walking';
    movement{2} = 'Straight Slow Running';
    movement{3} = 'Straight Fast Running';
    movement{4} = 'Curved Slow Running';
    movement{5} = 'V-Cut Movement';
    
    value = movement;
end
%======================================================================
%> @brief Get all names of the T-pose trials of OMC data
%>
%> @retval value is cell array with valid T-pose OMC trial names for each 
%> subject. It is of size numSubjects x 1. Each of this cell contains chars (
%> without file ending). 
%======================================================================
function value = getNPoseTrialNamesCell()
    trialNames = cell(11, 1);
    trialNames{1} = 'trial0002';
    trialNames{2} = 'trial0002';
    trialNames{3} = 'trial0002';
    trialNames{4} = 'trial0002';
    trialNames{5} = 'trial0002';
    trialNames{6} = 'trial0002';
    trialNames{7} = 'trial0002';
    trialNames{8} = 'trial0002';
    trialNames{9} = 'trial0002';
    trialNames{10}= 'trial0002';
    trialNames{11}= 'trial0002';
    value = trialNames;
end
%======================================================================
%> @brief Get all names of the trials of OMC data (used for cutting)
%>
%> @retval value is cell array with valid OMC trial names for each subject
%> for each movement. It is of size numSubjects x numMovements. Each of
%> this cell contains chars (without file ending). The names should match
%> to the order of getTriggerIndicesCell(). Currently, I am processing all
%> files for which the force plate(s) were hit correctly regardless of the
%> desired velocity.
%======================================================================
function value = getOMCTrialNamesCell()
    trialNames = cell(11, 5);
    % subject 2 
    trialNames{2, 1} = char({'trial0010'});  % straightwalking => Not used in the end
    trialNames{2, 2} = char({'trial0025', 'trial0026', 'trial0027'});  % straightslowrunning
    trialNames{2, 3} = char({'trial0042'});  % straightfastrunning => Not used in the end
    trialNames{2, 4} = char({'trial0100', 'trial0103', 'trial0105'});  % curvedslowrunning
    trialNames{2, 5} = char({'trial0122', 'trial0123', 'trial0130'});  % vcut
    % subject 3
    trialNames{3, 1} = char({});  % straightwalking => Not used in the end
    trialNames{3, 2} = char({'trial0028', 'trial0029', 'trial0036'});  % straightslowrunning
    trialNames{3, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{3, 4} = char({'trial0075', 'trial0077', 'trial0078'});  % curvedslowrunning
    trialNames{3, 5} = char({'trial0089', 'trial0095', 'trial0099'});  % vcut
    % subject 4
    trialNames{4, 1} = char({});  % straightwalking => Not used in the end
    trialNames{4, 2} = char({'trial0020', 'trial0030', 'trial0035'});  % straightslowrunning
    trialNames{4, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{4, 4} = char({'trial0063', 'trial0064', 'trial0067'});  % curvedslowrunning
    trialNames{4, 5} = char({'trial0079', 'trial0080', 'trial0081'});  % vcut
    % subject 5
    trialNames{5, 1} = char({});  % straightwalking => Not used in the end
    trialNames{5, 2} = char({'trial0031', 'trial0037', 'trial0043'});  % straightslowrunning
    trialNames{5, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{5, 4} = char({'trial0061', 'trial0071', 'trial0076'});  % curvedslowrunning
    trialNames{5, 5} = char({'trial0084', 'trial0085', 'trial0089'});  % vcut
    % subject 6
    trialNames{6, 1} = char({});  % straightwalking => Not used in the end
    trialNames{6, 2} = char({'trial0025', 'trial0029', 'trial0031'});  % straightslowrunning
    trialNames{6, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{6, 4} = char({'trial0076', 'trial0092', 'trial0093'});  % curvedslowrunning
    trialNames{6, 5} = char({'trial0103', 'trial0105', 'trial0107'});  % vcut
    % subject 7
    trialNames{7, 1} = char({});  % straightwalking => Not used in the end
    trialNames{7, 2} = char({'trial0026', 'trial0028', 'trial0031'});  % straightslowrunning
    trialNames{7, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{7, 4} = char({'trial0076', 'trial0077', 'trial0079'});  % curvedslowrunning
    trialNames{7, 5} = char({'trial0091', 'trial0101', 'trial0106'});  % vcut
    % subject 8
    trialNames{8, 1} = char({});  % straightwalking => Not used in the end
    trialNames{8, 2} = char({'trial0025', 'trial0026', 'trial0029'});  % straightslowrunning
    trialNames{8, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{8, 4} = char({'trial0061', 'trial0062', 'trial0064'});  % curvedslowrunning
    trialNames{8, 5} = char({'trial0076', 'trial0077', 'trial0091'});  % vcut
    % subject 9
    trialNames{9, 1} = char({});  % straightwalking => Not used in the end
    trialNames{9, 2} = char({'trial0030', 'trial0035', 'trial0036'});  % straightslowrunning
    trialNames{9, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{9, 4} = char({'trial0117', 'trial0119', 'trial0125'});  % curvedslowrunning
    trialNames{9, 5} = char({'trial0152', 'trial0155', 'trial0157'});  % vcut
    % subject 10
    trialNames{10, 1} = char({});  % straightwalking => Not used in the end
    trialNames{10, 2} = char({'trial0032', 'trial0036', 'trial0037'});  % straightslowrunning
    trialNames{10, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{10, 4} = char({'trial0079', 'trial0092', 'trial0127'});  % curvedslowrunning
    trialNames{10, 5} = char({'trial0104', 'trial0106', 'trial0110'});  % vcut
    % subject 11
    trialNames{11, 1} = char({});  % straightwalking => Not used in the end
    trialNames{11, 2} = char({'trial0026', 'trial0027', 'trial0028'});  % straightslowrunning
    trialNames{11, 3} = char({});  % straightfastrunning => Not used in the end
    trialNames{11, 4} = char({'trial0065', 'trial0067', 'trial0076'});  % curvedslowrunning
    trialNames{11, 5} = char({'trial0093', 'trial0094', 'trial0095'});  % vcut
    value = trialNames;
end