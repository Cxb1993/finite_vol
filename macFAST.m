function [C,fs,DELTA,STEPS] = macFAST(x,y,X,Y,dis,shocks,plot_her,Courant,TAU)
%Check arg number to see if we want to plot the schocks on results
    forloop = false;
    if nargin < 5
        shocks = 0; dis = 0;
        plot_her = true; Courant = 1;
    elseif nargin < 6
        shocks = 0; plot_her = true;Courant = 1;
    elseif nargin < 7
        plot_her = true;Courant = 1;
    elseif nargin == 9;
        forloop = true;
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
        newU = U; newU_ = U; TOL = 10^-10;
        t = 0; delta = 1; steps = 0;
        DELTA = zeros(5000,1);
        STEPS = zeros(5000,1);
    %set up grid for the animated plot (if so desire)
        if animate, plot_grid(x,y,X,Y,false);hold on; end;      
%------------------MAIN TIME STEPPING LOOP---------------------------------    
    %for tau=1:TAU%It want to manually set time steps
    while delta > TOL
            if mod(steps,200) == 0, 
                delta 
                steps
            end
            %if t>1.0*10^-2; break; end;
        %Set the slave cells
            U = set_slave(U);
        %Set E and F based on this adjusted U
            [E,F,rho,u,v,p] = recalc_EF(U,E,F);
        %Need a new dt each loop. 
            ua = abs(u);
            va = abs(v);
            c = (gamma*p./rho).^0.5;
            dt = min(min(Courant./((ua(2:JL-1,2:IL-1)/dx)+...
                (va(2:JL-1,2:IL-1)/dy_min)...
                +c(2:JL-1,2:IL-1)*(1/dx^2+1/dy_min^2)^0.5)));
            t = t+dt; steps = steps + 1;
        %Declare the dotted flux vector to be iterated with
            Ep_dot = zeros(size(U)); E_dot = Ep_dot;
            Fp_dot = Ep_dot; F_dot = Ep_dot;
        %Simplify indexing using simplifications a and b
            a = 2:JL-1; b = 2:IL-1;           
    %--------PREDICTOR-----------------------------------------------------
        %Calculate Dissapation Terms
            [dphi,dmhi,dphj,dmhj] = dis_terms2(u,v,c,p,N,dis); 
        %Begin the loop
             for i=1:4
             %Calculate dotted flux quantities
                Ep_dot(a,b,i) = E(a,3:IL,i).*N(a,b,1,2)... 
                              + F(a,3:IL,i).*N(a,b,2,2);
                E_dot(a,b,i)  = E(a,b,i).*N(a,b,1,4)...
                              + F(a,b,i).*N(a,b,2,4);
                Fp_dot(a,b,i) = E(3:JL,b,i).*N(a,b,1,3)...
                              + F(3:JL,b,i).*N(a,b,2,3);
                F_dot(a,b,i)  = E(a,b,i).*N(a,b,1,1)...
                              + F(a,b,i).*N(a,b,2,1);  
             %Calculate the predictor
                newU_(a,b,i) = U(a,b,i) - dt./V(a,b).*...
                  ((Ep_dot(a,b,i)-dphi.*(U(a,3:IL,i)-U(a,b,i))).*S(a,b,2)... 
                  +(E_dot(a,b,i)+dmhi.*(U(a,b,i)-U(a,1:IL-2,i))).*S(a,b,4)...
                  +(Fp_dot(a,b,i)-dphj.*(U(3:JL,b,i)-U(a,b,i))).*S(a,b,3)...
                  +(F_dot(a,b,i)+dmhj.*(U(a,b,i)-U(1:JL-2,b,i))).*S(a,b,1));
             end
        %Set the slave cells
            newU_ = set_slave(newU_);
        %Set E and F based on this adjusted U
            [E,F,~,u,v,p] = recalc_EF(newU_,E,F);
        %Declare dotted flux vectors to be used in Corrector
            Em_dot = zeros(size(U)); Fm_dot = Em_dot;            
    %--------Corrector-----------------------------------------------------
        %Calculate Dissapation Terms
            [dphi,dmhi,dphj,dmhj] = dis_terms4(u,v,c,p,N,dis); 
        %Begin the loop
            for i=1:4
            %Calculate dotted flux qunatities
                E_dot(a,b,i)  = E(a,b,i).*N(a,b,1,2)...
                              + F(a,b,i).*N(a,b,2,2);
                Em_dot(a,b,i) = E(a,1:IL-2,i).*N(a,b,1,4)...
                              + F(a,1:IL-2,i).*N(a,b,2,4);
                F_dot(a,b,i)  = E(a,b,i).*N(a,b,1,3)...
                              + F(a,b,i).*N(a,b,2,3);
                Fm_dot(a,b,i) = E(1:JL-2,b,i).*N(a,b,1,1)...
                              + F(1:JL-2,b,i).*N(a,b,2,1);   
           %Calculate the correctpr
               newU(a,b,i) = 0.5*(U(a,b,i) + newU_(a,b,i) - dt./V(a,b).*...
                  ((E_dot(a,b,i)-dphi.*(U(a,3:IL,i)-U(a,b,i))).*S(a,b,2)... 
                  +(Em_dot(a,b,i)+dmhi.*(U(a,b,i)-U(a,1:IL-2,i))).*S(a,b,4)...
                  +(F_dot(a,b,i)-dphj.*(U(3:JL,b,i)-U(a,b,i))).*S(a,b,3)...
                  +(Fm_dot(a,b,i)+dmhj.*(U(a,b,i)-U(1:JL-2,b,i))).*S(a,b,1)));
            end
        %Calculate delta 
                %disp(steps);
            if ~forloop
                delta = abs(U(2:JL-1,2:IL-1,4) - newU(2:JL-1,2:IL-1,4))...
                    ./U(2:JL-1,2:IL-1,4);
                delta = max(max(delta));
            elseif steps == TAU
                U = newU;
                break;
            end
        %Set U to the new values and we are ready to move on
            U = newU;
        %If the internal flag is set animate
            if animate,
               [C,fs]=plot_results(x,y,X,Y,U,shocks,true,true);
            end
        %DELTA STEPS
            DELTA(steps) = delta;
            STEPS(steps) = steps;
    end
    fprintf('Convergence took:\n\t %d seconds\n\t%d steps\n',t,steps);
%--------------------------------------------------------------------------
%POST-PROCESSING
%--------------------------------------------------------------------------
if ~animate
    [C,fs] = plot_results(x,y,X,Y,U,shocks,false,plot_her);
end
%--------------------------------------------------------------------------
%Auxillarily Functions
%--------------------------------------------------------------------------
    function [e] = energy(rho, p, u,v)
    gamma = 1.4;
    e = p/(gamma-1)+rho.*((u.^2+v.^2)/2);
    end
end