function [exact,shocks,fs] = exact_sol(x,y,X,Y,plot_her)
%check arg numbers
    fs = 0;
    if nargin < 5
        plot_her = true;
    end
%Declare Grid Size
    dim = size(x);
    grid_res = (dim(2)-2)/40;
    IL = 40*grid_res+2; 
    JL = 20*grid_res+2; 
    IS =  5*grid_res+1;
    dx = x(1,2)-x(1,1);
%conversion factor from degrees to rads
    d2r = pi/180;
%set flag for wheter of not we want to lpot the shocks on teh grid
    showshocks = true;
    if ~plot_her
        showshocks=false;
    end
%--------------------------------------------------------------------------
% FIND FLOW PARAMATERS IN EACH SECTION
%--------------------------------------------------------------------------
    %SECTION 1
        p1 = 10^5;
        rho1 = 1;
        M1 = 2.9;       
        theta = 10.94; 
        beta1_2 = tbM(theta,M1);
        M1n = sin(beta1_2*d2r)*M1;
    %SECTION 2
        [M2n, p2, rho2] = normal_shock(M1n,p1,rho1);
        M2 = M2n/(sin((beta1_2-theta)*d2r));
        %disp(M2); disp(p2); disp(rho2);
    %SECTION 3
        beta2_3 = tbM(theta,M2);
        M2n = sin(beta2_3*d2r)*M2;
        [M3n, p3, rho3] = normal_shock(M2n,p2,rho2);
        M3 = M3n/(sin((beta2_3-theta)*d2r));
        %disp(M3); disp(p3); disp(rho3);
    %SECTION 4
        beta3_4 = tbM(theta,M3);
        M3n = sin(beta3_4*d2r)*M3;
        [M4n, p4, rho4] = normal_shock(M3n,p3,rho3);
        M4 = M4n/(sin((beta3_4-theta)*d2r));
        %disp(M4); disp(p4); disp(rho4);
%--------------------------------------------------------------------------
% ASSIGN VALUES TO EACH POINT IN GRID
%--------------------------------------------------------------------------
    %Will accomplish this by makeing three lines of form y = mx+b
    %   where m is determined by beta and b is determined by known grid
    %   points
    %First create an equation to represent the top boundary
        mb = -tan(theta*d2r);
        bb = 1-mb*(dx*(IS-1));
    %FIRST DIVISION
        m1 = -tan(beta1_2*d2r);
        b1 = 1-m1*(dx*(IS-1));
        xcept = -b1/m1;        
    %SECOND DIVISION       
        m2 = tan((beta2_3-theta)*d2r);
        b2 = -m2*xcept;
    %THIRD DIISION
        xinter = (bb-b2)/(m2-mb);
        yinter = mb*xinter+bb;
        m3 = -tan(beta3_4*d2r);
        b3 = yinter-m3*xinter;
    %Optionally plot the shocks on the grid
        shock_xs = [dx*(IS-1),xcept,xinter, dx*(IL-1)];
        shock_ys = [1,0,yinter,m3*dx*(IL-1)+b3];
        shocks = [shock_xs;shock_ys];
        if showshocks 
            plotshocks(shock_xs,shock_ys);
        end
%Now create an array of appropriate dimension to hole the solutions    
        exact = zeros(dim(1),dim(2),3);
        S1 = [M1,p1,rho1]; S2 = [M2,p2,rho2]; 
        S3 = [M3,p3,rho3]; S4 = [M4,p4,rho4];
        C = zeros(dim(1)+1,dim(2)+1,3);
        for i=1:IL
            for j=1:JL
                %disp(strcat(num2str(i),'--',num2str(j)));
                test = [y(j,i)>=m1*x(j,i)+b1,y(j,i)>=m2*x(j,i)+b2...
                       ,y(j,i)>=m3*x(j,i)+b3];
                    if test(1)==false, exact(j,i,:)=S1;
                elseif test(2)== true, exact(j,i,:)=S2;
                elseif test(3)==false, exact(j,i,:)=S3;
                  else                 exact(j,i,:)=S4; 
                    end
                if i>1 && j>1 && i <IL && j <JL
                    C(j,i,:) = exact(j,i,:);
                end
            end
        end
%--------------------------------------------------------------------------
% PLOT RESULTS ON THE GRIDS
%--------------------------------------------------------------------------       
if plot_her
%Mach Number Plot
    fm = plot_grid(x,y,X,Y,false,false);hold on;
    title('Mach Number'); xlabel('x (m)');ylabel('y (m)');
    pcolor(X(2:JL,2:IL),Y(2:JL,2:IL),C(2:JL,2:IL,1)); shading FLAT;
    plotshocks(shock_xs,shock_ys);
    c=colorbar; A = ylabel(c,'Mach Number');  
    set(A,'Interpreter','Latex');
%Pressure Plot
    fp = plot_grid(x,y,X,Y,false,false);hold on;
    title('Pressure'); xlabel('x (m)');ylabel('y (m)');
    pcolor(X(2:JL,2:IL),Y(2:JL,2:IL),C(2:JL,2:IL,2)); shading FLAT;
    plotshocks(shock_xs,shock_ys);
    c=colorbar; A = ylabel(c,'Pressure (Pa)');
    set(A,'Interpreter','Latex');
%Density Plot
    fd = plot_grid(x,y,X,Y,false,false);hold on;
    title('Density'); xlabel('x (m)');ylabel('y (m)');
    pcolor(X(2:JL,2:IL),Y(2:JL,2:IL),C(2:JL,2:IL,3)); shading FLAT;
    plotshocks(shock_xs,shock_ys);
    c=colorbar; A =ylabel(c,'Density ($kg/m^3$)');
    set(A,'Interpreter','Latex');
fs = [fm fp fd];
end
%--------------------------------------------------------------------------
% AUX FUNCTIONS
%--------------------------------------------------------------------------
    function [] = plotshocks(xs,ys)
         hold on;
         plot([xs(1),xs(2)],[ys(1),ys(2)],'b','LineWidth',2);
         plot([xs(2),xs(3)],[ys(2),ys(3)],'b','LineWidth',2);
         plot([xs(3),xs(4)],[ys(3),ys(4)],'b','LineWidth',2.2);
    end
    
end