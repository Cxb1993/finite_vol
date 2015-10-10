function [U] = set_slave(U)
%Declare Grid Size
    dim = size(U);
    grid_res = (dim(2)-2)/40;
    IL = 40*grid_res+2; 
    JL = 20*grid_res+2; 
    IS =  5*grid_res+1;
%Set Slave Cells at exit
    U(:,IL,:) = 2*U(:,IL-1,:)-U(:,IL-2,:);
%Set Slave cells at top before shock
    U(JL,2:IS-1,1) = U(JL-1,2:IS-1,1); % density
    U(JL,2:IS-1,2) = U(JL-1,2:IS-1,2); % x velocity
    U(JL,2:IS-1,3) = -U(JL-1,2:IS-1,3);% y velocity
    U(JL,2:IS-1,4) = U(JL-1,2:IS-1,4); % energy
%Set slave cells at top after shock  
    %compute normal and tangental velocity want to set normal to negative
        theta = 10.94; 
        vn = U(JL-1,IS:IL,2)*sind(theta)+U(JL-1,IS:IL,3)*cosd(theta);
        vt = U(JL-1,IS:IL,2)*cosd(theta)-U(JL-1,IS:IL,3)*sind(theta);
    %Now set everyhting equal liek normal
        U(JL,IS:IL,2) = vt*cosd(theta)-vn*sind(theta);% x velocity
        U(JL,IS:IL,3) =-vt*sind(theta)-vn*cosd(theta);% y velocity
        U(JL,IS:IL,1) = U(JL-1,IS:IL,1); % density
        U(JL,IS:IL,4) = U(JL-1,IS:IL,4); % energy
%Set Slave cells at bottom
    U(1,:,1) = U(2,:,1); % density
    U(1,:,2) = U(2,:,2); % x velocity
    U(1,:,3) = -U(2,:,3);% y velocity
    U(1,:,4) = U(2,:,4); % energy            
end