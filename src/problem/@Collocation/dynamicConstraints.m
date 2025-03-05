%======================================================================
%> @file @Collocation/dynamicConstraints.m
%> @brief Collocation function to compute dynamic constraint
%> @details
%> Details: Collocation::dynamicConstraints()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding dynamic equilibrium
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = dynamicConstraints(obj,option,X)
%% check input parameter
if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
    error('Model states and controls and duration need to be stored in state vector X.')
end

%% state variable indices for semi-implicit Euler method
% https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
% Discretization is defined by ixSIE1 (state variables from the node before the time step)
% and ixSIE2 (state variables from the node after the time step).
% Two versions: SIE A and SIE B. Both conserve energy.
% SIE A appeared easier to solve in tests.  This makes sense because it uses backward
% Euler for muscle equations. SIE B code is commented out and provided for comparison.
% SIE A: 
ixSIE2 = sort([obj.model.extractState('qdot'); obj.model.extractState('s'); obj.model.extractState('a')]);
ixSIE1 = setdiff(1:obj.model.nStates, ixSIE2);
% SIE B
% ixSIE1 = sort([obj.model.extractState('q')]);
% ixSIE2 = setdiff(1:obj.model.nStates, ixSIE1);


%% compute demanded output
nNodesDur = obj.nNodesDur;
h = X(obj.idx.dur)/(nNodesDur-1);
nconstraintspernode = obj.model.nConstraints;

if strcmp(option,'confun')
    output = zeros(nconstraintspernode*(nNodesDur-1),1);
    
    % dynamic equations must be zero
    for iNode=1:(nNodesDur-1)
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        x1 = X(obj.idx.states(:,iNode));
        x2 = X(obj.idx.states(:,iNode+1));
        xd =(x2-x1)/h;
        u1 = X(obj.idx.controls(:,iNode));
        u2 = X(obj.idx.controls(:,iNode+1));
        
        if strcmp(obj.Euler,'BE')
            output(ic) = obj.model.getDynamics(x2,xd,u2);	% backward Euler discretization
        elseif strcmp(obj.Euler,'ME')
            % we're using u2 instead of (u1+u2)/2 because u is
            % open loop and it makes no difference except u2
            % will converge more easily
            output(ic) = obj.model.getDynamics((x1+x2)/2,xd,u2);
        elseif strcmp(obj.Euler,'SIE')
            % semi-implicit Euler is energy neutral
			% take specified states ixSIE2 from x2, all other states from x1
			x1(ixSIE2) = x2(ixSIE2);
			output(ic) = obj.model.getDynamics(x1,xd,u2);
        % higher-order methods
		elseif strcmp(obj.Euler, 'HS1') %First order Hermite simpson
			f1 = obj.model.getDynamics(x1, xd, u1)+xd;
			f2 = obj.model.getDynamics(x2, xd, u2)+xd;
			
			xc = 1/2*(x1+x2)+h/8*(f1-f2);
			uc = 1/2*(u1+u2); %Is this described in the paper?
			fc = obj.model.getDynamics(xc, xd, uc) + xd;
			output(ic) = xd-1/6*(f1+4*fc+f2);
        elseif strcmp(obj.Euler, 'HS2') %Second order Hermite simpson
            nq = 1:obj.model.nDofs; %this should be done nicer using the state table
            nv = obj.model.nDofs+1:obj.model.nDofs*2;
            nmus = obj.model.nDofs*2+(1:obj.model.nMus*2);
            nc = obj.model.nDofs*2+obj.model.nMus*2+(1:obj.model.nCPs*4);
            
            q1 = x1(nq);
            v1 = x1(nv);
            mus1 = x1(nmus);
            c1 = x1(nc);
            v2 = x2(nv);
            mus2 = x2(nmus);
            c2 = x2(nc);

            f1 = obj.model.getDynamics(x1, xd, u1)+xd;
			f2 = obj.model.getDynamics(x2, xd, u2)+xd;
			
            g1 = f1(nv);
            g2 = f2(nv);

            qc = q1+h/32*(13*v1+3*v2)+h^2/192*(11*g1-5*g2);
            vc = 1/2*(v1+v2)+h/8*(g1-g2);
            musc= 1/2*(mus1+mus2)+h/8*(f1(nmus)-f2(nmus));
            cc= 1/2*(c1+c2)+h/8*(f1(nc)-f2(nc));

            xc = [qc;vc;musc;cc];
            uc = 1/2*(u1+u2); %Is this described in the paper?
            fc = obj.model.getDynamics(xc,xd,uc)+xd;
            gc = fc(nv);

            output(ic(nq)) = xd(nq) - v1 - h/6*(g1+2*gc);
            output(ic(nv)) = xd(nv) - 1/6*(g1+4*gc+g2);
            output(ic(nmus)) = xd(nmus) - 1/6* (f1(nmus)+4*fc(nmus)+f2(nmus));
            output(ic(nc)) = xd(nc) - 1/6* (f1(nc)+4*fc(nc)+f2(nc));
        end
    end
elseif strcmp(option,'jacobian')
    output = spalloc(nconstraintspernode*(nNodesDur-1),length(X),obj.Jnnz);
    
    for iNode = 1:(nNodesDur-1)
        ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
        
        ix1 = obj.idx.states(:,iNode);
        ix2 = obj.idx.states(:,iNode+1);
        iu1 = obj.idx.controls(:,iNode);
        iu2 = obj.idx.controls(:,iNode+1);
        x1 = X(ix1);
        x2 = X(ix2);
        xd = (x2-x1)/h;
        u1 = X(iu1);
        u2 = X(iu2);
        
        if strcmp(obj.Euler,'BE')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x2,xd,u2);
            output(ic,ix1) = -dfdxdot'/h;
            output(ic,ix2) = dfdx' + dfdxdot'/h;
            output(ic,iu2) = dfdu';
            % derivative of constraints with respect to duration (because h is duration/(N-1))
            output(ic,obj.idx.dur) = -dfdxdot' * (x2-x1) / h^2 / (nNodesDur-1);
        elseif strcmp(obj.Euler,'ME')
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics((x1+x2)/2,xd,u2);
            output(ic,ix1) = dfdx'/2 - dfdxdot'/h;
            output(ic,ix2) = dfdx'/2 + dfdxdot'/h;
            output(ic,iu2) = dfdu';
            % derivative of constraints with respect to duration (because h is duration/(N-1))
            output(ic,obj.idx.dur) = -dfdxdot' * (x2-x1) / h^2 / (nNodesDur-1);
        elseif strcmp(obj.Euler,'SIE')
			x1(ixSIE2) = x2(ixSIE2);
            [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x1,xd,u2);
            output(ic,ix1) = -dfdxdot'/h;
            output(ic,ix2) = dfdxdot'/h;
            output(ic,ix1(ixSIE1)) = output(ic,ix1(ixSIE1)) + dfdx(ixSIE1,:)';
            output(ic,ix2(ixSIE2)) = output(ic,ix2(ixSIE2)) + dfdx(ixSIE2,:)';
            output(ic,iu2) = dfdu';
            % derivative of constraints with respect to duration (because h is duration/(N-1))
            output(ic,obj.idx.dur) = -dfdxdot' * (x2-x1) / h^2 / (nNodesDur-1);
        % higher-order methods
		elseif strcmp(obj.Euler, 'HS1') %First order Hermite simpson
			[notf1, df1dx, df1dxdot, df1du1] = obj.model.getDynamics(x1,xd,u1);
            [notf2, df2dx, df2dxdot, df2du2] = obj.model.getDynamics(x2,xd,u2);
            f1 = notf1 + xd;
            df1dx1 = df1dx' - df1dxdot'/h - 1/h*speye(length(x1));
            df1dx2 = df1dxdot'/h + 1/h*speye(length(x2));
            df1dh = df1dxdot'*-(x2-x1)/h^2-(x2-x1)/h^2;

            f2 = notf2 + xd;
            df2dx1 = -df2dxdot'/h - 1/h*speye(length(x1));
            df2dx2 = df2dx' + df2dxdot'/h + 1/h*speye(length(x2)) ;
            df2dh = df2dxdot'*-(x2-x1)/h^2-(x2-x1)/h^2;

            xc = 1/2*(x1+x2)+h/8*(f1-f2);
            dxcdx1 = 1/2*speye(length(x1)) + h/8*(df1dx1 - df2dx1);
            dxcdx2 = 1/2*speye(length(x2)) + h/8*(df1dx2 - df2dx2);
            dxcdu1 = h/8*(df1du1');
            dxcdu2 = h/8*(-df2du2');
            dxcdh = 1/8*(f1-f2) + h/8*(df1dh - df2dh);

			uc = 1/2*(u1+u2); %Is this described in the paper?

			[~, dfcdxc, dfcdxdot, dfcduc] = obj.model.getDynamics(xc,xd,uc);
            % fc = notfc + xd;           
            dfcdx1 = dfcdxc'*dxcdx1 - dfcdxdot'/h - 1/h*speye(length(x1)) ;
            dfcdx2 = dfcdxc'*dxcdx2 + dfcdxdot'/h + 1/h*speye(length(x2));
            dfcdu1 = 1/2*dfcduc' + dfcdxc'*dxcdu1;
            dfcdu2 = 1/2*dfcduc' + dfcdxc'*dxcdu2;
            dfcdh = dfcdxc'*dxcdh + dfcdxdot'*-(x2-x1)/h^2 - (x2-x1)/h^2;

			% output(ic) = xd-1/6*(f1+4*fc+f2);
            output(ic,ix1) = -1/h*speye(length(x1))-1/6*(df1dx1 + 4*dfcdx1 + df2dx1);
            output(ic,ix2) = 1/h*speye(length(x2))-1/6*(df1dx2 + 4*dfcdx2 + df2dx2);
            output(ic,iu1) = -1/6*(df1du1' + 4*dfcdu1);
            output(ic,iu2) = -1/6*(4*dfcdu2 + df2du2');

            % derivative of constraints with respect to duration (because h is duration/(N-1))
            output(ic,obj.idx.dur) = 1/(nNodesDur-1)*(-(x2-x1)/h^2 - 1/6*(df1dh + 4*dfcdh + df2dh));
        elseif strcmp(obj.Euler, 'HS2') %Second order Hermite simpson
            nq = 1:obj.model.nDofs; %this should be done nicer using the state table
            nv = obj.model.nDofs+1:obj.model.nDofs*2;
            nmus = obj.model.nDofs*2+(1:obj.model.nMus*2);
            nc = obj.model.nDofs*2+obj.model.nMus*2+(1:obj.model.nCPs*4);
            
            q1 = x1(nq);
            v1 = x1(nv);
            mus1 = x1(nmus);
            c1 = x1(nc);
            v2 = x2(nv);
            mus2 = x2(nmus);
            c2 = x2(nc);

            [notf1, df1dx, df1dxdot, df1du1] = obj.model.getDynamics(x1,xd,u1);
            [notf2, df2dx, df2dxdot, df2du2] = obj.model.getDynamics(x2,xd,u2);
            f1 = notf1 + xd;
            df1dx1 = df1dx' - df1dxdot'/h - 1/h*speye(length(x1));
            df1dx2 = df1dxdot'/h + 1/h*speye(length(x2));
            df1dh = df1dxdot'*-(x2-x1)/h^2-(x2-x1)/h^2;

            f2 = notf2 + xd;
            df2dx1 = -df2dxdot'/h - 1/h*speye(length(x1));
            df2dx2 = df2dx' + df2dxdot'/h + 1/h*speye(length(x2)) ;
            df2dh = df2dxdot'*-(x2-x1)/h^2-(x2-x1)/h^2;

            g1 = f1(nv); 
            dg1dx1 = df1dx1(nv,:);
            dg1dx2 = df1dx2(nv,:);
            dg1du1 = df1du1(:,nv)';
            dg1dh = df1dh(nv,:);
            
            g2 = f2(nv);
            dg2dx1 = df2dx1(nv,:);
            dg2dx2 = df2dx2(nv,:);
            dg2du2 = df2du2(:,nv)';
            dg2dh = df2dh(nv,:);

            qc = q1+h/32*(13*v1+3*v2)+h^2/192*(11*g1-5*g2);
            dqcdq1 = speye(obj.model.nDofs) + h^2/192*(11*dg1dx1(:,nq)-5*dg2dx1(:,nq));
            dqcdq2 = h^2/192*(11*dg1dx2(:,nq)-5*dg2dx2(:,nq));
            dqcdv1 = speye(obj.model.nDofs)*13*h/32 + h^2/192*(11*dg1dx1(:,nv)-5*dg2dx1(:,nv));
            dqcdv2 =  speye(obj.model.nDofs)*3*h/32 + h^2/192*(11*dg1dx2(:,nv)-5*dg2dx2(:,nv));
            dqcdmus1 = h^2/192*(11*dg1dx1(:,nmus)-5*dg2dx1(:,nmus));
            dqcdmus2 = h^2/192*(11*dg1dx2(:,nmus)-5*dg2dx2(:,nmus));
            dqcdc1 = h^2/192*(11*dg1dx1(:,nc)-5*dg2dx1(:,nc));
            dqcdc2 = h^2/192*(11*dg1dx2(:,nc)-5*dg2dx2(:,nc));
            dqcdx1 = [dqcdq1 dqcdv1 dqcdmus1 dqcdc1];
            dqcdx2 = [dqcdq2 dqcdv2 dqcdmus2 dqcdc2];
            dqcdu1 = h^2/192*(11*dg1du1);
            dqcdu2 = h^2/192*(-5*dg2du2);
            dqcdh = 1/32*(13*v1+3*v2)+2*h/192*(11*g1-5*g2) + h^2/192*(11*dg1dh-5*dg2dh);

            vc = 1/2*(v1+v2)+h/8*(g1-g2);
            dvcdq1 = h/8*(dg1dx1(:,nq)-dg2dx1(:,nq));
            dvcdq2 = h/8*(dg1dx2(:,nq)-dg2dx2(:,nq));
            dvcdv1 = speye(obj.model.nDofs)/2 + h/8*(dg1dx1(:,nv)-dg2dx1(:,nv));
            dvcdv2 = speye(obj.model.nDofs)/2 + h/8*(dg1dx2(:,nv)-dg2dx2(:,nv));
            dvcdmus1 = h/8*(dg1dx1(:,nmus)-dg2dx1(:,nmus));
            dvcdmus2 = h/8*(dg1dx2(:,nmus)-dg2dx2(:,nmus));
            dvcdc1 = h/8*(dg1dx1(:,nc)-dg2dx1(:,nc));
            dvcdc2 = h/8*(dg1dx2(:,nc)-dg2dx2(:,nc));
            dvcdx1 = [dvcdq1 dvcdv1 dvcdmus1 dvcdc1];
            dvcdx2 = [dvcdq2 dvcdv2 dvcdmus2 dvcdc2];
            dvcdu1 = h/8*dg1du1;
            dvcdu2 = -h/8*dg2du2;
            dvcdh = 1/8*(g1-g2) +  h/8*(dg1dh-dg2dh);

            musc= 1/2*(mus1+mus2)+h/8*(f1(nmus)-f2(nmus));
            dmuscdq1 = h/8*(df1dx1(nmus,nq)-df2dx1(nmus,nq));
            dmuscdq2 = h/8*(df1dx2(nmus,nq)-df2dx2(nmus,nq));
            dmuscdv1 = h/8*(df1dx1(nmus,nv)-df2dx1(nmus,nv));
            dmuscdv2 = h/8*(df1dx2(nmus,nv)-df2dx2(nmus,nv));
            dmuscdmus1 = speye(obj.model.nMus*2)/2 + h/8*(df1dx1(nmus,nmus)-df2dx1(nmus,nmus));
            dmuscdmus2 = speye(obj.model.nMus*2)/2 + h/8*(df1dx2(nmus,nmus)-df2dx2(nmus,nmus));
            dmuscdc1 = h/8*(df1dx1(nmus,nc)-df2dx1(nmus,nc));
            dmuscdc2 = h/8*(df1dx2(nmus,nc)-df2dx2(nmus,nc));
            dmuscdx1 = [dmuscdq1 dmuscdv1 dmuscdmus1 dmuscdc1];
            dmuscdx2 = [dmuscdq2 dmuscdv2 dmuscdmus2 dmuscdc2];
            dmuscdu1 = h/8*df1du1(:,nmus)';
            dmuscdu2 = -h/8*df2du2(:,nmus)';
            dmuscdh = 1/8*(f1(nmus)-f2(nmus)) +  h/8*(df1dh(nmus,:)-df2dh(nmus,:));

            cc= 1/2*(c1+c2)+h/8*(f1(nc)-f2(nc));
            dccdq1 = h/8*(df1dx1(nc,nq)-df2dx1(nc,nq));
            dccdq2 = h/8*(df1dx2(nc,nq)-df2dx2(nc,nq));
            dccdv1 = h/8*(df1dx1(nc,nv)-df2dx1(nc,nv));
            dccdv2 = h/8*(df1dx2(nc,nv)-df2dx2(nc,nv));
            dccdmus1 = h/8*(df1dx1(nc,nmus)-df2dx1(nc,nmus));
            dccdmus2 = h/8*(df1dx2(nc,nmus)-df2dx2(nc,nmus));
            dccdc1 = h/8*(df1dx1(nc,nc)-df2dx1(nc,nc))+ speye(obj.model.nCPs*4)/2;
            dccdc2 = h/8*(df1dx2(nc,nc)-df2dx2(nc,nc))+ speye(obj.model.nCPs*4)/2;
            dccdx1 = [dccdq1 dccdv1 dccdmus1 dccdc1];
            dccdx2 = [dccdq2 dccdv2 dccdmus2 dccdc2];
            dccdu1 = h/8*df1du1(:,nc)';
            dccdu2 = -h/8*df2du2(:,nc)';
            dccdh = 1/8*(f1(nc)-f2(nc)) +  h/8*(df1dh(nc,:)-df2dh(nc,:));

            xc = [qc;vc;musc;cc];
            dxcdx1 = [dqcdx1;dvcdx1;dmuscdx1;dccdx1];
            dxcdx2 = [dqcdx2;dvcdx2;dmuscdx2;dccdx2];
            dxcdu1 = [dqcdu1;dvcdu1;dmuscdu1;dccdu1];
            dxcdu2 = [dqcdu2;dvcdu2;dmuscdu2;dccdu2];
            dxcdh = [dqcdh;dvcdh;dmuscdh;dccdh];

            uc = 1/2*(u1+u2); %Is this described in the paper?
            
            [notfc, dfcdxc, dfcdxdot, dfcduc] = obj.model.getDynamics(xc,xd,uc);
            fc = notfc+xd;
            dfcdx1 = dfcdxc'*dxcdx1 - dfcdxdot'/h - 1/h*speye(length(x1));
            dfcdx2 = dfcdxc'*dxcdx2 + dfcdxdot'/h + 1/h*speye(length(x2));
            dfcdu1 = 1/2*dfcduc' + dfcdxc'*dxcdu1;
            dfcdu2 = 1/2*dfcduc' + dfcdxc'*dxcdu2;
            dfcdh = dfcdxc'*dxcdh + dfcdxdot'*-(x2-x1)/h^2 - (x2-x1)/h^2;

            gc = fc(nv);
            dgcdx1 = dfcdx1(nv,:);
            dgcdx2 = dfcdx2(nv,:);
            dgcdu1 = dfcdu1(nv,:);
            dgcdu2 = dfcdu2(nv,:);
            dgcdh = dfcdh(nv,:);
            
            %output(ic(nq)) = xd(nq) - v1 - h/6*(g1+2*gc);
            output(ic(nq),ix1(nq)) = -1/h*speye(length(nq));
            output(ic(nq),ix1(nv)) = -speye(length(nv));
            output(ic(nq),ix1) = output(ic(nq),ix1) - h/6*(dg1dx1 + 2*dgcdx1);

            output(ic(nq),ix2(nq)) = 1/h*speye(length(nq));
            output(ic(nq),ix2) = output(ic(nq),ix2) - h/6*(dg1dx2 + 2*dgcdx2);

            output(ic(nq),iu1) = -h/6*(dg1du1 + 2*dgcdu1);
            output(ic(nq),iu2) = -h/6*(2*dgcdu2);

            % derivative with respect to duration (because h is duration/(N-1))
            %Do the first part for all indices for speed
            output(ic,obj.idx.dur) = 1/(nNodesDur-1)*-(x2-x1)/h^2;
            output(ic(nq),obj.idx.dur) = output(ic(nq),obj.idx.dur) + 1/(nNodesDur-1)*(- 1/6*(g1+2*gc) - h/6*(dg1dh+2*dgcdh));

            %output(ic(nv)) = xd(nv) - 1/6*(g1+4*gc+g2);
            output(ic(nv),ix1(nv)) = -1/h*speye(length(nq));
            output(ic(nv),ix2(nv)) = 1/h*speye(length(nq));
            output(ic(nv),ix1) = output(ic(nv),ix1) - 1/6*(dg1dx1 + 4*dgcdx1 + dg2dx1);
            output(ic(nv),ix2) = output(ic(nv),ix2) - 1/6*(dg1dx2 + 4*dgcdx2 + dg2dx2);
            output(ic(nv),iu1) = -1/6*(dg1du1 + 4*dgcdu1);
            output(ic(nv),iu2) = -1/6*(4*dgcdu2 + dg2du2);

            % derivative of constraints with respect to duration (because h is duration/(N-1))
            output(ic(nv),obj.idx.dur) = output(ic(nv),obj.idx.dur) + 1/(nNodesDur-1)*(- 1/6*(dg1dh + 4*dgcdh + dg2dh));

            %output(ic(nmus)) = xd(nmus) - 1/6* (f1(nmus)+4*fc(nmus)+f2(nmus));
            output(ic(nmus),ix1(nmus)) = -1/h*speye(length(nmus));
            output(ic(nmus),ix2(nmus)) = 1/h*speye(length(nmus));
            output(ic(nmus),ix1) = output(ic(nmus),ix1) - 1/6*(df1dx1(nmus,:) + 4*dfcdx1(nmus,:) + df2dx1(nmus,:));
            output(ic(nmus),ix2) = output(ic(nmus),ix2) - 1/6*(df1dx2(nmus,:) + 4*dfcdx2(nmus,:) + df2dx2(nmus,:));
            output(ic(nmus),iu1) = -1/6*(df1du1(:,nmus)' + 4*dfcdu1(nmus,:));
            output(ic(nmus),iu2) = -1/6*(4*dfcdu2(nmus,:) + df2du2(:,nmus)');

            % derivative of constraints with respect to duration (because h is duration/(N-1))
            output(ic(nmus),obj.idx.dur) = output(ic(nmus),obj.idx.dur) + 1/(nNodesDur-1)*(- 1/6*(df1dh(nmus,:) + 4*dfcdh(nmus,:) + df2dh(nmus,:)));
            
            %output(ic(nc)) = xd(nc) - 1/6* (f1(nc)+4*fc(nc)+f2(nc));
            output(ic(nc),ix1(nc)) = -1/h*speye(length(nc));
            output(ic(nc),ix2(nc)) = 1/h*speye(length(nc));
            output(ic(nc),ix1) = output(ic(nc),ix1) - 1/6*(df1dx1(nc,:) + 4*dfcdx1(nc,:) + df2dx1(nc,:));
            output(ic(nc),ix2) = output(ic(nc),ix2) - 1/6*(df1dx2(nc,:) + 4*dfcdx2(nc,:) + df2dx2(nc,:));
            output(ic(nc),iu1) = -1/6*(df1du1(:,nc)' + 4*dfcdu1(nc,:));
            output(ic(nc),iu2) = -1/6*(4*dfcdu2(nc,:) + df2du2(:,nc)');

            % derivative of constraints with respect to duration (because h is duration/(N-1))
            output(ic(nc),obj.idx.dur) = output(ic(nc),obj.idx.dur) + 1/(nNodesDur-1)*(- 1/6*(df1dh(nc,:) + 4*dfcdh(nc,:) + df2dh(nc,:)));

        end
        
    end
else
    error('Unknown option.');
end
end
