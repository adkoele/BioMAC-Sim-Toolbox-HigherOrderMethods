%======================================================================
%> @file @Collocation/dynamicConstraints.m
%> @brief Collocation function to compute dynamic constraints using higher
%>          order methods
%> @details
%> Details: Collocation::dynamicConstraintsHO()
%>
%> @author Alexander Weiss
%> @date March, 2025
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding dynamic equilibrium
%>
%> @param obj           Collocation class object
%>
% @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = dynamicConstraintsRK(obj,option,X, discMeth)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'xk1_states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('Model states, intermediate states and controls and duration need to be stored in state vector X.')
end

%% compute demanded output
nNodesDur = obj.nNodesDur;
h = X(obj.idx.dur)/(nNodesDur-1);

%% get Butcher Tableau
[A, b, c] = DiscMeth.getButcherTableau(discMeth);

nStages = size(A,1); % Number of RK stages
n_k = size(A,2);

nconstraintspernode = obj.model.nConstraints*(1+nStages);

if strcmp(option,'confun')
    output = zeros((nconstraintspernode)*(nNodesDur-1),1);
    % dynamic equations must be zero
    for iNode=1:(nNodesDur-1)
        ix1 = obj.idx.states(:,iNode);
        ix2 = obj.idx.states(:,iNode+1);

        iu1 = obj.idx.controls(:,iNode);
        iu2 = obj.idx.controls(:,iNode+1);

        u1 = X(iu1);
        u2 = X(iu2);

        x1 = X(ix1);
        x2 = X(ix2);
        xd = (x2-x1)/h;
        nx = size(x1,1);

        ic = (1:obj.model.nConstraints) + (iNode-1)*nconstraintspernode;

        for i_k = 1:nStages
            ic_xk{i_k} = ((i_k-1)*nconstraintspernode/(nStages+1)+size(ic,2)+1:nconstraintspernode*((i_k)/(nStages+1))+size(ic,2)) + (iNode-1)*(nconstraintspernode); %indices of constraints of iNode in c
            xk_Name = sprintf('xk%d_states', i_k);
            ixkx{i_k} = obj.idx.(xk_Name)(:,iNode);
            kxnStat(:,i_k) = X(ixkx{i_k});
        end

        xdStates = cell(1,nStages);

        for i = 1:nStages
            xdStates{i} = (kxnStat(:,i) - x1) / h;
        end
        for i = 1:nStages
            weightedSum = zeros(nx, 1); 
            for k = 1:n_k
                weightedSum = weightedSum + A(i, k) * xdStates{k};
            end
            xkStates{i} = x1 + h * weightedSum;               
            output(ic_xk{i}) = obj.model.getDynamics(xkStates{i}, xdStates{i}, u1*(1-c(i)) + c(i)*u2);
        end

        sum_xd = zeros(size(x1, 1),1);
        for i_k = 1:nStages
            k_Name = sprintf('xk%d_states', i_k);
            ixk{i_k} = obj.idx.(k_Name)(:,iNode);
            % Variante
            %knStat(:,i_k) = X(ixk{i_k});
            sum_xd = sum_xd + b(i_k) * X(ixk{i_k});%obj.model.getDynamics(X(ixk{i_k}),xd,u2);
        end    
        output(ic) = x2 - sum_xd;

    end
elseif strcmp(option,'jacobian')
    output = spalloc((nconstraintspernode)*(nNodesDur-1),length(X),obj.Jnnz);
    for iNode = 1:(nNodesDur-1)
        ix1 = obj.idx.states(:,iNode);
        ix2 = obj.idx.states(:,iNode+1);

        iu1 = obj.idx.controls(:,iNode);
        iu2 = obj.idx.controls(:,iNode+1);

        x1 = X(ix1);
        x2 = X(ix2);
        xd = (x2-x1)/h;

        u1 = X(iu1);
        u2 = X(iu2);

        
        ic = (1:obj.model.nConstraints) + (iNode-1)*nconstraintspernode;

        for i_k = 1:nStages
            ic_xk{i_k} = ((i_k-1)*nconstraintspernode/(nStages+1)+size(ic,2)+1:nconstraintspernode*((i_k)/(nStages+1))+size(ic,2)) + (iNode-1)*(nconstraintspernode); %indices of constraints of iNode in c
            xk_Name = sprintf('xk%d_states', i_k);
            ixkx{i_k} = obj.idx.(xk_Name)(:,iNode);
            kxnStat(:,i_k) = X(ixkx{i_k});
        end

        nx = size(x1,1);
 
        % Store xk-states dynamically
        kStates = cell(1, n_k);
        for k = 1:n_k
            xkStates{k} = kxnStat(:,k);
        end

        % Initialize matrices
        dx = cell(nStages, n_k);  
        xdStates = cell(1,nStages);

        for i = 1:nStages
            xdStates{i} = (kxnStat(:,i) - x1) / h;
            for k = 1:n_k      
                dx{i,k} = + h * A(i,k) * eye(nx) / h;
            end
        end

        % Derivatives of xi w.r.t x1 and itself (assume x1 is constant)
        dxidx1 = repmat({eye(nx)}, 1, nStages);  % dx2dx1, dx3dx1, ...
        dxidx2 = repmat({zeros(nx)}, 1, nStages);

        for i = 1:nStages
            weightedSum = zeros(nx, 1);  % initialize with same size as x1
            for k = 1:n_k               
                weightedSum = weightedSum + A(i, k) * xdStates{k};
            end           
            xkStates{i} = x1 + h * weightedSum;
                           
            [~, dfdx{i}, dfdxdot{i}, dfdu{i}] = obj.model.getDynamics(xkStates{i}, xdStates{i}, u1*(1-c(i)) + c(i)*u2); 
  
            % Derivatives w.r.t. each k-stateStates{i}
            for k = 1:n_k
                dkcx{1,i}.dk{i,k} = dfdx{i}' * dx{i,k};  
            end
            dkcx{1,i}.dk{i,i} = + dkcx{1,i}.dk{i,i} + dfdxdot{i}' / h;
            % Derivative w.r.t. h
            weightedSum_k = zeros(nx, 1);
            weightedSum_x1 = zeros(nx, nx);

            for k = 1:n_k
                weightedSum_x1 = weightedSum_x1 + A(i,k) * - eye(nx) / h;
                output(ic_xk{i}, ixkx{k}) = dkcx{1,i}.dk{i,k};
            end
            
            dkcx{i}.dh =  - dfdxdot{i}' * (xdStates{i}) / (h) / (nNodesDur-1);
            output(ic_xk{i}, ix1) = dfdx{i}' * (eye(nx) + h * weightedSum_x1) - dfdxdot{i}' / h;        

            output(ic_xk{i}, iu1) = dfdu{i}'*(1-c(i));
            output(ic_xk{i}, iu2) = dfdu{i}'*c(i);
            output(ic_xk{i}, obj.idx.dur) = dkcx{i}.dh;
        end

        for i_k = 1:nStages
            k_Name = sprintf('xk%d_states', i_k);
            ixk{i_k} = obj.idx.(k_Name)(:,iNode);
            output(ic,ixk{i_k}) = - b(i_k) * eye(size(x2,1));
        end
        output(ic,ix2) = eye(size(x2,1));
    end
else
    error('Unknown option.');
end
end
