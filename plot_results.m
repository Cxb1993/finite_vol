function [C,fs] = plot_results(x,y,X,Y,U,shocks,animate,plot_her)
    fs = 0;
%Declare Grid Size
    dim = size(U);
    grid_res = (dim(2)-2)/40;
    IL = 40*grid_res+2; 
    JL = 20*grid_res+2; 
    IS =  5*grid_res+1;
%declare gamma and if we are showing shocks
    gamma = 1.4;
    showshocks = true;
    if shocks(1) == 0
        showshocks = false;
    end
%Extract M p and rho from U
    RHO = U(:,:,1);
    u = U(:,:,2)./RHO;
    v = U(:,:,3)./RHO;
    PRESSURE = pressure(RHO,U(:,:,4),u,v);
    c = (gamma*PRESSURE./RHO).^0.5;
    M = ((u.^2 + v.^2).^0.5)./c;
%Set up color plot
    C = zeros(dim(1)+1,dim(2)+1,3);
    for i=2:IL-1
        for j=2:JL-1
                C(j,i,:) = [M(j,i), PRESSURE(j,i),RHO(j,i)];
        end
    end
    if ~plot_her, return; end
if ~animate
    %----------------------------------------------------------------------
    % PLOT RESULTS ON THE GRIDS
    %----------------------------------------------------------------------
    %Mach Number Plot
        fm = plot_grid(x,y,X,Y,false,false);hold on;
        title('Mach Number'); xlabel('x (m)');ylabel('y (m)');
        pcolor(X(2:JL,2:IL),Y(2:JL,2:IL),C(2:JL,2:IL,1)); shading FLAT;
        if showshocks, plotshocks(shocks(1,:),shocks(2,:)); end
        c=colorbar; A = ylabel(c,'Mach Number');  
        set(A,'Interpreter','Latex');
    %Pressure Plot
        fp = plot_grid(x,y,X,Y,false,false);hold on;
        title('Pressure'); xlabel('x (m)');ylabel('y (m)');
        pcolor(X(2:JL,2:IL),Y(2:JL,2:IL),C(2:JL,2:IL,2)); shading FLAT;
        if showshocks, plotshocks(shocks(1,:),shocks(2,:)); end
        c=colorbar; A = ylabel(c,'Pressure (Pa)');
        set(A,'Interpreter','Latex');
    %Density Plot
        fd = plot_grid(x,y,X,Y,false);hold on;
        title('Density'); xlabel('x (m)');ylabel('y (m)');
        pcolor(X(2:JL,2:IL),Y(2:JL,2:IL),C(2:JL,2:IL,3)); shading FLAT;
        if showshocks, plotshocks(shocks(1,:),shocks(2,:)); end
        c=colorbar; A =ylabel(c,'Density $kg/m^3$');
        set(A,'Interpreter','Latex');
    fs = [fm, fp, fd];
else
    %Density Plot
        title('Density'); xlabel('x (m)');ylabel('y (m)');
        pcolor(X(2:JL,2:IL),Y(2:JL,2:IL),C(2:JL,2:IL,3)); shading FLAT;
        if showshocks, plotshocks(shocks(1,:),shocks(2,:)); end
        c=colorbar; A =ylabel(c,'Density $kg/m^3$');
        set(A,'Interpreter','Latex');
        drawnow();
    fs = 0;
end
%--------------------------------------------------------------------------
%Auxillarily Functions
%--------------------------------------------------------------------------
    function [p] = pressure(rho, e, u,v)
    gamma = 1.4;
    p = (gamma-1)*(e-rho.*(u.^2+v.^2)/2);
    end
    function [] = plotshocks(xs,ys)
         hold on;
         plot([xs(1),xs(2)],[ys(1),ys(2)],'b','LineWidth',2);
         plot([xs(2),xs(3)],[ys(2),ys(3)],'b','LineWidth',2);
         plot([xs(3),xs(4)],[ys(3),ys(4)],'b','LineWidth',2.2);
    end
end