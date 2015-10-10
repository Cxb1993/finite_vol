function [U] = macdis(x,y,X,Y,dis,shocks,TAU)
%Check arg number to see if we want to plot the schocks on results
    plot_her = true;
    if nargin < 5
        shocks = 0; dis = 0;
        TAU = 500;
    elseif nargin < 6
        shocks = 0; TAU = 500;
    elseif nargin < 7
        TAU = 500;
    else
        plot_her = false;
    end    
%OPTION: Animate plot with this flag
    animate = false;
%--------------------------------------------------------------------------
%PRE-PROCESSING: Settig up structures for numerical method
%--------------------------------------------------------------------------
%Declare Grid Size
    dim = size(x);
    grid_res = (dim(2)-2)/40;
    IL = 40*grid_res+2; 
    JL = 20*grid_res+2; 
    IS =  5*grid_res+1;
    dx = x(1,2)-x(1,1);
    dy_min = Y(2,IL)-Y(1,IL);
%Givens. ICs and BCs
    p1 = 10^5;
    rho1 = 1;
    M1 = 2.9;        
    gamma = 1.4;
%Some more important needed quantites for numerical implemetation. 
    c1 = sqrt(gamma*p1/rho1);
    u1 = M1*c1;
    v1 = 0;
    e1 = energy(rho1,p1,u1,v1);
    U1 = [rho1;rho1*u1;rho1*v1;e1];
    E1 = [rho1*u1;rho1*u1^2+p1;rho1*u1*v1;(e1+p1)*u1];
    F1 = [rho1*v1;rho1*u1*v1;rho1*v1^2+p1;(e1+p1)*v1];
%Make structures to hold each cells pertinent values.
    %geometry based - Area(V)m side lengths(S), normals(N)
        N = zeros(dim(1),dim(2),2,4);
        S = zeros(dim(1),dim(2),4);
        V = zeros(dim);
    %flow based - flux vectors E and F. Also U
        E = zeros(dim(1),dim(2),4);
        F = E;
        U = F;
    %Set up grid for these quantities
        for i=1:dim(2)
            for j=1:dim(1)
                [n,s,v] = cell_normals(X,Y,[j,i]);
                N(j,i,:,:) = n;
                S(j,i,:) = s;
                V(j,i) = v;
                E(j,i,:) = E1;
                F(j,i,:) = F1;
                U(j,i,:) = U1;
            end
        end
%--------------------------------------------------------------------------
%NUMERICAL METHOD: MACCORMACK
%--------------------------------------------------------------------------
    %Some variables needed for iteration
        newU = U; newU_ = U;
        t = 0; delta = 1; steps = 0;
    %set up grid for the animated plot (if so desire)
        if animate, plot_grid(x,y,X,Y,false,false);hold on; end;
    %Loop it    
    for tau=1:TAU%It want to manually set time steps
    %while delta >0.0001
        %Set the slave cells
            U = set_slave(U);
        %Set E and F based on this adjusted U
            [E,F,rho,u,v,p] = recalc_EF(U,E,F);
        %Need a new dt each loop. 
            u = abs(u);
            v = abs(v);
            c = (gamma*p./rho).^0.5;
            dt = min(min(1./((u(2:JL-1,2:IL-1)/dx)+...
                (v(2:JL-1,2:IL-1)/dy_min)...
                +c(2:JL-1,2:IL-1)*(1/dx^2+1/dy_min^2)^0.5)));
            t = t+dt; steps = steps + 1;
        %Loop over each cell oce to make predictor
            for i=2:IL-1
            for j=2:JL-1
            %Calculate the flux quantities
                Ep_dot = E(j,i+1,:)*N(j,i,1,2)+F(j,i+1,:)*N(j,i,2,2);
                E_dot  = E(j,i,:)*N(j,i,1,4)+F(j,i,:)*N(j,i,2,4);
                Fp_dot = E(j+1,i,:)*N(j,i,1,3)+F(j+1,i,:)*N(j,i,2,3);
                F_dot = E(j,i,:)*N(j,i,1,1)+F(j,i,:)*N(j,i,2,1);
            %Calculate Dissapation Terms - not sure what to do with uprime
                [dphi,dmhi,dphj,dmhj] = dis_terms(u,v,c,p,[j,i],dis);              
            %Calculate the predictor
                newU_(j,i,:) = U(j,i,:) - dt/V(j,i)*...
                    (S(j,i,2)*(Ep_dot-dphi*(U(j,i+1,:)-U(j,i,:)))...
                   + S(j,i,4)*(E_dot +dmhi*(U(j,i,:)-U(j,i-1,:)))...
                   + S(j,i,3)*(Fp_dot-dphj*(U(j+1,i,:)-U(j,i,:)))...
                   + S(j,i,1)*(F_dot +dmhj*(U(j,i,:)-U(j-1,i,:))));        
            end    
            end
        %Set the slave cells
            newU_ = set_slave(newU_);
        %Set E and F based on this adjusted U
            [E,F,~,u,v,p] = recalc_EF(newU_,E,F);
        %Loop over each cell to do corrector
            for i=2:IL-1
            for j=2:JL-1
            %Calculate the flux quantities
                E_dot =E(j,i,:)*N(j,i,1,2)+F(j,i,:)*N(j,i,2,2);
                Em_dot=E(j,i-1,:)*N(j,i,1,4)+F(j,i-1,:)*N(j,i,2,4);
                F_dot =E(j,i,:)*N(j,i,1,3)+F(j,i,:)*N(j,i,2,3);
                Fm_dot=E(j-1,i,:)*N(j,i,1,1)+F(j-1,i,:)*N(j,i,2,1);
            %Calculate Dissapation Terms - not sure what to do with uprime
                [dphi,dmhi,dphj,dmhj] = dis_terms(u,v,c,p,[j,i],dis); 
            %Calculate the corrector
                newU(j,i,:) = 0.5*(U(j,i,:)+newU_(j,i,:)-dt/V(j,i)*...  
                    (S(j,i,2)*(E_dot -dphi*(U(j,i+1,:)-U(j,i,:)))...
                   + S(j,i,4)*(Em_dot+dmhi*(U(j,i,:)-U(j,i-1,:)))...
                   + S(j,i,3)*(F_dot -dphj*(U(j+1,i,:)-U(j,i,:)))...
                   + S(j,i,1)*(Fm_dot+dmhj*(U(j,i,:)-U(j-1,i,:))))); 
            end
            end   
        %Set U to the new values and we are ready to move on
            delta = mean([abs(U(ceil(14*JL/22),ceil(25*IL/42),1)...
                            -newU(ceil(14*JL/22),ceil(25*IL/42),1)),
                         abs(U(JL-1,IL-1,1)-newU(JL-1,IL-1,1))]);
            U = newU;
            if animate,
               plot_results(x,y,X,Y,U,shocks,true);
            end
    end
    fprintf('Convergence took:\n\t %d seconds\n\t%d steps\n',t,steps);
%--------------------------------------------------------------------------
%POST-PROCESSING
%--------------------------------------------------------------------------
if ~animate && plot_her
    plot_results(x,y,X,Y,U,shocks,false);
end
%--------------------------------------------------------------------------
%Auxillarily Functions
%--------------------------------------------------------------------------
    function [e] = energy(rho, p, u,v)
    gamma = 1.4;
    e = p/(gamma-1)+rho.*((u.^2+v.^2)/2);
    end
end