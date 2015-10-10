function [bulk_error,fs] = compareToExact(Cexact,C,x,plot_her)
%check arg numbers
    if nargin < 4
        plot_her = true;
    end
%Declare Grid Size
    dim = size(x);
    grid_res = (dim(2)-2)/40;
    IL = 40*grid_res+2; 
    JL = 20*grid_res+2;
%Plot the characteristic lines Top Middle and end
        bulk_error = 0;
    %Determine plotting domain
        span = 1:IL-2;
        top = JL-2;
        bot = 1;
        mid = round(JL/2);
    %top
    if plot_her
        ft=figure; hold on;
        title('Exact and Numerical Solutions for the Top Elements');
        xlabel('x(m)'); ylabel('Mach Number');
        plot(x(top,span),Cexact(top+1,span+1,1));
        plot(x(top,span),C(top+1,span+1,1),'.r','MarkerSize',15);
    end
%         bulk_error = bulk_error+...
%             sum(abs(Cexact(top+1,span+1,1)-C(top+1,span+1,1)));
    %mid
    if plot_her
        fm=figure; hold on;
        title('Exact and Numerical Solutions for the Middle Elements');
        xlabel('x(m)'); ylabel('Mach Number');
        plot(x(mid,span),Cexact(mid+1,span+1,1));
        plot(x(mid,span),C(mid+1,span+1,1),'.r','MarkerSize',15);
    end
%         bulk_error = bulk_error+...
%             sum(abs(Cexact(mid+1,span+1,1)-C(mid+1,span+1,1)));
    %bot
    if plot_her
        fb=figure; hold on;
        title('Exact and Numerical Solutions for the Bottom Elements');
        xlabel('x(m)'); ylabel('Mach Number');
        plot(x(bot,span),Cexact(bot+1,span+1,1));
        plot(x(bot,span),C(bot+1,span+1,1),'.r','MarkerSize',15);
    end
%         bulk_error = bulk_error+...
%             sum(abs(Cexact(bot+1,span+1,1)-C(bot+1,span+1,1)));
        
    bulk_error = bulk_error + ...
        sum(sum(abs(Cexact(bot+1:top+1,span+1)-C(bot+1:top+1,span+1))));
    if plot_her, fs = [ft,fm,fb]; 
    else fs = 0;
    end
end