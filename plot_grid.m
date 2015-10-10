function [f] = plot_grid(x,y,x2,y2,show_inode,show_mesh)
%INPUT ARG CHECK
switch nargin
    case 5
        show_mesh = true;
    case 4
        show_inode = true; show_mesh = true;
    otherwise
end
%REFIND SOME NEEDED VALUES 
    %Declare Grid Size
        dim = size(x);
        grid_res = (dim(2)-2)/40;
        IL = 40*grid_res+2; 
        JL = 20*grid_res+2; 
        IS =  5*grid_res+1; 
%ATUALLY PLOT THE THING
    %Start a new figure and plot the 
        f = figure; hold on;
        title('Numerical Grid');
        xlabel('x (m)');
        ylabel('y (m)');
    %plot the boundaries with xs
        bnode=plot(x(2:dim(1)-1,1),y(2:dim(1)-1,1),'kx');
        plot(x(2:dim(1)-1,dim(2)),y(2:dim(1)-1,dim(2)),'kx');
        plot(x(1,:),y(1,:),'kx');
        plot(x(dim(1),:),y(dim(1),:),'kx')
    %plot the internal nodes with .s
    if show_inode
        inode=plot(x(2:dim(1)-1,2:dim(2)-1),y(2:dim(1)-1,2:dim(2)-1),'k.');
        legend([inode(1), bnode],'Internal Nodes','Boundary Nodes');
    end
    %plot the separation lines
        %Vertical
            for i=1:IL+1
                if show_mesh
                    plot([x2(1,i),x2(dim(1)+1,i)],...
                         [y2(1,i),y2(dim(1)+1,i)],'k');
                elseif i<=2 || i>= IL
                    plot([x2(1,i),x2(dim(1)+1,i)],...
                         [y2(1,i),y2(dim(1)+1,i)],'k');
                end
            end
        %Horizontal (ish)
            for j=1:JL+1
                if show_mesh
                    plot([x2(j,1),x2(j,IS),x2(j,IL+1)],...
                         [y2(j,1),y2(j,IS),y2(j,IL+1)],'k');
                elseif j<=2 || j>= JL
                    plot([x2(j,1),x2(j,IS),x2(j,IL+1)],...
                         [y2(j,1),y2(j,IS),y2(j,IL+1)],'k');
                end
            end
end