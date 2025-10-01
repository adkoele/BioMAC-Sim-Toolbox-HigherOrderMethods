%======================================================================
%> @file @Collocation/dynamicConstraintsRK.m
%> @brief Collocation function to compute dynamic constraints using higher
%>          order methods
%> @details
%> Details: Collocation::dynamicConstraintsRK()
%>
%> @author Alexander Weiss
%> @date March, 2025
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding dynamic equilibrium
%>
%> @param obj           Collocation class object
%>
%> @param option         String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%> @param discMeth       String: Name of discretization method (e.g. BE or RIIa-3)
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
n_k = size(A,2);     % Second dimension of matrix, typically equal to number of RK stages

nconstraintspernode = obj.model.nConstraints*(1+nStages);

if strcmp(option,'confun')
    output = zeros((nconstraintspernode)*(nNodesDur-1),1);

    nx = size(obj.idx.states,1);  % Assuming consistent size across all nodes
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
        sum_xd = zeros(size(x1, 1),1);

        ic = (1:obj.model.nConstraints) + (iNode-1)*nconstraintspernode;

        for i_k = 1:nStages
            % Indices of RK stages
            ic_xk(i_k,:) = ((i_k-1)*nconstraintspernode/(nStages+1)+size(ic,2)+1:nconstraintspernode*((i_k)/(nStages+1))+size(ic,2)) + (iNode-1)*(nconstraintspernode); %indices of constraints of iNode in c
            % RK stages names
            xk_Name = sprintf('xk%d_states', i_k);
            ixkx(i_k,:) = obj.idx.(xk_Name)(:,iNode);
            % Position of RK stages in X
            kxnStat(:,i_k) = X(ixkx(i_k,:));
            % Numerical differntiation of RK stages
            xdStates(i_k,:) = (kxnStat(:,i_k) - x1) / h;
        end

        for i = 1:nStages
            weightedSum = zeros(nx, 1);
            for k = 1:n_k
                weightedSum = weightedSum + A(i, k) * xdStates(k,:)';
            end
            % Definition of RK stages according to Butcher Tableaus
            xkStates(:,i) = x1 + h * weightedSum;
            % Linear interpolation of u, depending on c of Butcher Tableau
            % for specific method
            output(ic_xk(i,:)) = obj.model.getDynamics(xkStates(:,i), xdStates(i,:)', u1*(1-c(i)) + c(i)*u2);

            % Weighted sum depending on Butcher Tableau
            sum_xd = sum_xd + b(i) * kxnStat(:,i);
        end
        output(ic) = x2 - sum_xd;
    end

  
elseif strcmp(option,'jacobian')

    output = spalloc((nconstraintspernode)*(nNodesDur-1),length(X),obj.Jnnz);
    for iNode = 1:(nNodesDur-1)
    
        % Indices
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
        nx = size(x1,1);
        sum_xd = zeros(size(x1, 1),1);

        for i_k = 1:nStages
            ic_xk(i_k,:) = ((i_k-1)*nconstraintspernode/(nStages+1)+size(ic,2)+1:nconstraintspernode*((i_k)/(nStages+1))+size(ic,2)) + (iNode-1)*(nconstraintspernode); %indices of constraints of iNode in c
            xk_Name = sprintf('xk%d_states', i_k);
            ixkx(:,i_k) = obj.idx.(xk_Name)(:,iNode);
            kxnStat(:,i_k) = X(ixkx(:,i_k));
            xdStates(:,i_k) = (kxnStat(:,i_k) - x1) / h;
        end

        for i = 1:nStages
            weightedSum = zeros(nx, 1);  % initialize with same size as x1
            for k = 1:n_k
                weightedSum = weightedSum + A(i, k) * xdStates(:,k);
            end
            xkStates(:,i) = x1 + h * weightedSum;

            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(xkStates(:,i), xdStates(:,i), u1*(1-c(i)) + c(i)*u2);

            weightedSum_x1 = zeros(nx, nx);
            % Derivatives w.r.t. each RK state
            for k = 1:n_k
                dx = + h * A(i,k) * eye(nx) / h;
                dkcxdk = dfdx' * dx;
                weightedSum_x1 = weightedSum_x1 + A(i,k) * - eye(nx) / h;
                output(ic_xk(i,:), ixkx(:,k)) = dkcxdk;
            end
            output(ic_xk(i,:), ixkx(:,i)) = output(ic_xk(i,:), ixkx(:,i)) + dfdxdot' / h;
            % Derivative w.r.t. h
            dkcxdh =  - dfdxdot' * (xdStates(:,i)) / (h) / (nNodesDur-1);
            output(ic_xk(i,:), ix1) = dfdx' * (eye(nx) + h * weightedSum_x1) - dfdxdot' / h;
            output(ic_xk(i,:), iu1) = dfdu'*(1-c(i));
            output(ic_xk(i,:), iu2) = dfdu'*c(i);
            output(ic_xk(i,:), obj.idx.dur) = dkcxdh;
            output(ic,ixkx(:,i)) = - b(i) * eye(size(x2,1));
        end
        output(ic,ix2) = eye(size(x2,1));
    end
else
    error('Unknown option.');
end
end
