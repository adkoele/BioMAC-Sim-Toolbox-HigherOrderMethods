% ======================================================================
%> @file Collocation.m
%> @brief Matlab class for defining an optimal control of musculoskeletal
%> model
%
%> @author Eva Dorschky, Marlies Nitschke
%> @date November, 2017
% ======================================================================

% ======================================================================
%> @brief The class for defining an optimal control of musculoskeletal
% ======================================================================
classdef Collocation_RK < Collocation

    properties
        
    end

    methods
        %======================================================================
        %> @brief Constructor setting default Collocation options.
        %>
        %> @param   model     Model class object
        %> @param   nNodes    (optional) Double: Number of collocation nodes
        %> @param   Euler     (optional) String: Discretization method: backward euler 'BE' or midpoint euler 'ME' or semi-implicit euler 'SIE
        %> @param   logfile   (optional) String: Logfile including the path. The input can be skipped with an empty string.
        %> @param   plotLog   (optional) Bool: If true, plot log in call of Collocation.objfun() (stick figure, objective, constrains) (default: 0)
        %>
        %> @retval    obj   Collocation class object
        %======================================================================
        function [obj] = Collocation_RK(model,nNodes,Euler,logfile,plotLog)
            obj = obj@Collocation(model,nNodes,'BE',logfile,plotLog);%varargin{:});
            tic
            
            if nargin < 2
                obj.nNodes = 50;
            elseif nNodes <= 0
                error('The number of collocation nodes is invalid.');
            elseif nNodes > 1600
                answer = input('The number of collocation nodes is very high and the simulation will likely be slow. Are you sure you want to continue? Press y/n',"s");
                if strcmp(answer,'y')
                    obj.nNodes = nNodes;
                else
                    return
                end
            else
                obj.nNodes = nNodes;
            end

            if nargin < 3
                obj.Euler = 'BE';
            elseif strcmp(Euler,'BE') || strcmp(Euler,'ME') || strcmp(Euler,'SIE') || strcmp(Euler,'LIIIc-2') || strcmp(Euler,'RIIa-2') || strcmp(Euler,'LIIIc-3') || strcmp(Euler,'RIIa-3')
                obj.Euler = Euler;
            else
                error('The discretization method is not in the valid range.');
            end
        end

      
        %======================================================================
        %> @brief Function to define an initial guess
        %>
        %> @details
        %> The initial guess is initialized by adding the state variables. This function can be used to
        %> overwrite the initial guess, e.g. by loading an initial guess of a previous
        %> result. The input parameter 'init' can be used to define how the
        %> initial guess should be computed. If init = 'mid', the midpoint
        %> between the lower and upper bounds of the state vector X are
        %> used. If init = 'random', X is intialized by drawing random
        %> values from a uniform distribution between the lower and upper
        %> bound. If init contains the path to a result file, the result
        %> file is loaded. If neccessary, the result is resampled to the
        %> size of the current collocation problem.
        %>
        %> @code
        %> problem.makeinitialguess('mid');
        %> @endcode
        %>
        %> @param   obj         Collocation class object
        %> @param   init        String: Parsing the demanded initialguess
        %>                      OR
        %>                      Double vector: X which fits exactly the size of
        %>                      the current problem (no tests and no resampling)
        %======================================================================
        function makeinitialguess(obj,init)

            if isnumeric(init)
                % init is X vector
                obj.initialguess.X = init;
                obj.initialguess.info.type = 'X';

            else
                % init is string
                switch init
                    case 'mid'
                        obj.initialguess.X = (obj.X_lb + obj.X_ub) / 2; % use mid point between lower and upper bound
                    case 'random'
                        rng('shuffle'); % shuffle the generator: https://de.mathworks.com/help/matlab/math/why-do-random-numbers-repeat-after-startup.html
                        obj.initialguess.X = obj.X_lb + (obj.X_ub-obj.X_lb).*rand(size(obj.X_lb));  % use random point from uniform distribution between lower and upper bound
                    case 'random_smart'
                        load('result');
                        state_names = fieldnames(obj.idx);
                        for iname = 1:length(state_names)
                            name = state_names{iname};
                            if ~isfield(result.problem.idx,name)
                                warning(['The variable ''',name,''' does not exist in initial guess file. The initial guess of this variable could not updated.']);
                                continue
                            elseif size(obj.idx.(name),1) ~= size(result.X(result.problem.idx.(name)),1);
                                warning(['The variable ''',name,''' has a different numbers of states. The initial guess of this variable could not updated.']);
                                continue
                            end
                            %> @todo What if symmetric and not symmetric
                            %> movement is loaded?
                            if size(obj.idx.(name),2) ~= size(result.X(result.problem.idx.(name)),2);
                                obj.initialguess.X(obj.idx.(name)) = resampleX(result.X(result.problem.idx.(name)),size(obj.idx.(name),2)); % function resampleX can be overloaded and adpated by the user
                            else
                                obj.initialguess.X(obj.idx.(name)) = result.X(result.problem.idx.(name));
                            end
                            obj.initialguess.info = result.info; %save information of result which is used as initial guess, e.g. for warmstart
                        end
                        rand_nums = 0.5 + (1.5-0.5).*rand(size(obj.initialguess.X));
                        obj.initialguess.X = obj.initialguess.X.*rand_nums;
                        % reinitialize mex for current model, as mex function was
                        % initialized with result model by loading result calling load(init,'result');
                        %
                        %> @todo: current hack without changing Model Class, is there
                        %> a better way to make sure that always the right model is
                        %> initialized?
                        if ~eq(obj.model,result.problem.model)
                            obj.model.saveobj; % saveobj sets obj.init = 0
                            obj.model = obj.model.loadobj(obj.model); % loadobj calls initMex when obj.init = 0
                        end
                    otherwise
                        try
                            load(init,'result');
                        catch
                            error('Make initialguess failed. Check whether result file exists and is in the search path.');
                        end
                        state_names = fieldnames(obj.idx);
                        for iname = 1:length(state_names)
                            name = state_names{iname};
                            if ~isfield(result.problem.idx,name)
                                if contains(name,'xk')
                                    if size(obj.idx.(name),2) ~= size(result.X(result.problem.idx.states),2);
                                        obj.initialguess.X(obj.idx.(name)) = resampleX(result.X(result.problem.idx.states),size(obj.idx.(name),2)); % function resampleX can be overloaded and adpated by the user
                                    else
                                        obj.initialguess.X(obj.idx.(name)) = result.X(result.problem.idx.states);
                                    end
                                else
                                    warning(['The variable ''',name,''' does not exist in initial guess file. The initial guess of this variable could not updated.']);
                                end
                                continue
                            elseif size(obj.idx.(name),1) ~= size(result.X(result.problem.idx.(name)),1);
                                warning(['The variable ''',name,''' has a different numbers of states. The initial guess of this variable could not updated.']);
                                continue
                            end
                            %> @todo What if symmetric and not symmetric
                            %> movement is loaded?
                            if size(obj.idx.(name),2) ~= size(result.X(result.problem.idx.(name)),2)
                                obj.initialguess.X(obj.idx.(name)) = resampleX(result.X(result.problem.idx.(name)),size(obj.idx.(name),2)); % function resampleX can be overloaded and adpated by the user
                            else
                                obj.initialguess.X(obj.idx.(name)) = result.X(result.problem.idx.(name));
                            end
                            obj.initialguess.info = result.info; %save information of result which is used as initial guess, e.g. for warmstart
                        end
                        % reinitialize mex for current model, as mex function was
                        % initialized with result model by loading result calling load(init,'result');
                        %
                        %> @todo: current hack without changing Model Class, is there
                        %> a better way to make sure that always the right model is
                        %> initialized?
                        if ~eq(obj.model,result.problem.model)
                            obj.model.saveobj; % saveobj sets obj.init = 0
                            obj.model = obj.model.loadobj(obj.model); % loadobj calls initMex when obj.init = 0
                        end
                end
                obj.initialguess.info.type = init;
            end
        end    
    end
end
%======================================================================
%> @brief Function for resampling state vector for makeinitialguess
%>
%> @param   x         Matrix with initial guess
%> @param   N         Number of collocation nodes
%> @retval  xrep      Matrix with resampled initial guess
%======================================================================
function xrep = resampleX(x,N)
if size(x,2) == 1
    xrep = repmat(x,1,N);
else
    xq = (0:N-1)/N * (size(x,2)-1) + 1;
    xrep = interp1(x',xq)' ;
end
end


