function [x, y, x2, y2] = gen_grid(grid_res)
%ACTUAL NODAL POSITIONS
    %Declare Grid Size
        IL = 40*grid_res+2; 
        JL = 20*grid_res+2; 
        IS =  5*grid_res+1; 
    %Declare Geometry of Problem
        H = 1; 
        L = 3.2;
        theta = 10.94;
    %Determine dx adn dy
        dx = L/(IL-IS);
        dy = zeros(1,IL);
        dh = tan(theta*pi()/180)*dx;
        dy(1:(IS-1)) = H/(JL-2);
        for i=IS:length(dy);
            dy(i) = (H-dh*(0.5+i-IS))/(JL-2);
        end
    %set up x and y coordinates in a matrix 
        x = meshgrid(0.5*dx:dx:IL*dx,1:JL);
        y = x;
        for i=1:IL
            y_change_vector = dy(i)*(-1:1:JL-2)';
            y(:,i) = (0.5*dy(i)+dy(1))*ones(JL,1)+y_change_vector;
        end
    %Shift plot down so the horixontal line corresponds with y=0
        y = y - ones(size(x))*dy(1);
%GRID CONTAINING THE BOARDERS FORE EACH OF TEH ABOVE NODES
    %set up boundary grid
        x2 = meshgrid(0:dx:(IL)*dx,1:JL+1);
        y2 = x2;
        dy = zeros(1,IL+1);
        dy(1:IS) = H/(JL-2);
        y2_2nd = dy(1);
        for i=IS+1:length(dy);
            dy(i) = (H-dh*(i-IS))/(JL-2);
        end
        for i=1:IL+1
            y2_change_vector = dy(i)*(-1:1:JL-1)';
            y2(:,i) = y2_2nd*ones(JL+1,1)+y2_change_vector;
        end
    %Shift plot down so the horixontal line corresponds with y=0
        y2 = y2 - ones(size(x2))*dy(1);
end