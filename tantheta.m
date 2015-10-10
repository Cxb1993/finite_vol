function [tant] = tantheta(M,beta)
gamma = 1.4;
tant = 2*(1/tan(beta))*(M^2*(sin(beta))^2-1)...
       /(M^2*(gamma+cos(2*beta))+2);
end