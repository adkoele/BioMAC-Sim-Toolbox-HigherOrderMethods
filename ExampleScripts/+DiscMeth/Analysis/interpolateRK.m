%======================================================================
%> @file    +DiscMeth/interpolateRK.m
%> @brief   Interpolate states, controls, and residuals between collocation nodes.
%>
%> @details
%> This function performs interpolation of simulation results obtained with 
%> Runge-Kutta based collocation methods. Given the discrete values of state, 
%> state derivatives, control signals, and dynamic residuals, the function 
%> generates smoothly interpolated trajectories using a specified number of 
%> interpolation points per interval.  
%>
%> @param   resultRunning   Struct: Simulation results structure containing 
%>                                   states, state derivatives, controls, and 
%>                                   dynamic residuals from collocation.
%> @param   numInterp       Int: Number of interpolation points to be 
%>                                   generated between collocation nodes.
%>
%> @retval  x_all_interp        Double array: Interpolated state trajectories.
%> @retval  xdot_all_interp     Double array: Interpolated state derivatives.
%> @retval  residual_all_interp Double array: Interpolated dynamic residuals.
%> @retval  u_all_interp        Double array: Interpolated control inputs.
%>
%>
%> @author  Alexander Weiss
%> @date    August, 2025
%======================================================================



function [x_all_interp, xdot_all_interp, residual_all_interp, u_all_interp] = interpolateRK(resultRunning, numInterp)
    X = resultRunning.X;
    nNodes = resultRunning.problem.nNodesDur;
    h = X(resultRunning.problem.idx.dur)/(nNodes-1);
    
    for iNode = 1:resultRunning.problem.nNodes-1 
        x1(:,iNode) = X(resultRunning.problem.idx.states(:,iNode))';
        u1(:,iNode) = X(resultRunning.problem.idx.controls(:,iNode))';
        u2(:,iNode) = X(resultRunning.problem.idx.controls(:,iNode+1))'; 
        nx = size(x1,1);
    
    
        [A,b,c] = DiscMeth.getButcherTableau(resultRunning.problem.Euler);
    
        nNodesDur = resultRunning.problem.nNodesDur;
        h = X(resultRunning.problem.idx.dur)/(nNodesDur-1);
    
        nStages = size(A,1); % Number of RK stages
        n_k = size(A,2);
     
        for i_k = 1:nStages
            xk_Name = sprintf('xk%d_states', i_k);
            ixkx{i_k} = resultRunning.problem.idx.(xk_Name)(:,iNode);
            kxnStat(:,i_k) = X(ixkx{i_k});
            xkn(:,iNode,i_k) = kxnStat(:,i_k);
        end
    
        for i = 1:nStages
            xdStates(:,i) = (xkn(:,iNode, i) - x1(:,iNode)) / h;
            fk(iNode,:, i) = resultRunning.problem.model.getDynamics(xkn(:,iNode, i), xdStates(:,i), u2(:,iNode));
            vk(iNode,:, i) = xdStates(:,i); 
        end
    end
    
    nx = size(x1,1);
    nInterp = numInterp; 
    nTotalInterp = (resultRunning.problem.nNodes ) * nInterp; 
    
    x_all_interp = zeros(nx, nTotalInterp);         % interpolated states
    xdot_all_interp = zeros(nx, nTotalInterp);      % interpolated derivatives
    residual_all_interp = zeros(nx, nTotalInterp);  % interpolation defects
    u_all_interp = zeros(size(u1,1), nTotalInterp); % interpolated controls
    
    
    for iNode = 1:resultRunning.problem.nNodes-1 %!!!!
        nInterp = numInterp;
        tau_interp = linspace(0, 1-(1/nInterp), nInterp);
        nx = size(x1,1);
    
        % Extract stage derivatives: k1, k2
        k1 = vk(iNode,:,1)'; % size nx x 1
        if nStages > 1
            k2 = vk(iNode,:,2)';
            if nStages > 2
                k3 = vk(iNode,:,3)';
            end
        end
    
        offset = (iNode - 1) * nInterp;
    
        for i_tau = 1:nInterp
            tau = tau_interp(i_tau) * h;
    
            % Interpolation polynomial (Radau IIA specific)
            x_tau = x1(:,iNode) + k1 * tau;
            xdot_tau = k1;
    
            if nStages == 2
                x_tau = x_tau + (k2 - k1)*tau^2 / (2*h);
                xdot_tau = xdot_tau + (k2-k1) * tau / h;
    
                if contains(resultRunning.problem.Euler, 'RIIa')
                    x_tau = x1(:,iNode) + 0.5*(3*k1 - k2) * tau + 0.75 * (k2 - k1) * tau^2 / h ;
                    xdot_tau = 0.5*(3*k1 - k2) + 1.5*(k2 - k1) * tau / h;
                end
            end
    
            if nStages == 3
                x_tau = x_tau - (3*k1 - 4*k2 + k3)*tau^2 / (2*h) + 2*(k1 - 2*k2 + k3)*tau^3 / (3*h^2);
                xdot_tau = xdot_tau - (3*k1 - 4*k2 + k3) * tau / (h) + 2*(k1 - 2*k2 + k3)*tau^2 / (h^2);
                if contains(resultRunning.problem.Euler, 'RIIa')
                    b = (3*sqrt(6)  * (k1-k2) + 2*(k1+k2+k3))/6;
                    c = (13*sqrt(6) * (k2-k1) + 8*(k1+k2-2*k3))/(12*h);
                    d = (5*sqrt(6)  * (k1-k2) - 5*(k1+k2-2*k3))/(9*h^2);
                    x_tau = x1(:,iNode) + b * tau + c * tau^2 + d * tau^3;
                    xdot_tau = b + 2 * c * tau + 3 * d * tau^2;
                end
            end
    
            % Linearly interpolate control (or use other interpolation if needed)
            u_tau = (1 - tau) * u1(:,iNode) + tau * u2(:,iNode);
    
            % Compute residual
            f_tau = resultRunning.problem.model.getDynamics(x_tau, xdot_tau, u_tau);
            residual = f_tau;
    
            % Store each in full matrix
            idx = offset + i_tau;
            x_all_interp(:, idx) = x_tau;
            xdot_all_interp(:, idx) = xdot_tau;
            residual_all_interp(:, idx) = residual;
            u_all_interp(:, idx) = u_tau;
        end
    end
end