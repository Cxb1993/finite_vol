function [normals,areas,volume] = cell_normals(X,Y,pt)
%get the cell we are intrested in
    j = pt(1); i = pt(2);
    cell_X = [X(j,i),X(j,i+1),X(j+1,i+1),X(j+1,i)];
    cell_Y = [Y(j,i),Y(j,i+1),Y(j+1,i+1),Y(j+1,i)];
%calculate the normals and areas (2D-Length)
    O = [1,2,3,4,1];
    normals = zeros(2,4);
    areas = zeros(1,4);
    for j=1:4;
        dx = cell_X(O(j+1)) - cell_X(O(j));
        dy = cell_Y(O(j+1)) - cell_Y(O(j));
        areas(j) = sqrt(dx^2+dy^2);
        n = [dy;-dx]/areas(j);
        normals(:,j) = n;        
    end
%calculate the volume of the cell (2D-Area)
    S1 = det([cell_X(1), cell_Y(1), 1 ;
              cell_X(3), cell_Y(3), 1 ;
              cell_X(4), cell_Y(4), 1 ]);
    S2 = det([cell_X(1), cell_Y(1), 1 ;
              cell_X(2), cell_Y(2), 1 ;
              cell_X(3), cell_Y(3), 1 ]);
    volume = 0.5*(abs(S1)+abs(S2));
end