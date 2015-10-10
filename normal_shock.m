function [M2n, p2, rho2] = normal_shock(M1n,p1, rho1)
    gamma = 1.4;
    M2n = sqrt((M1n^2*(gamma-1)+2)/(2*gamma*M1n^2-(gamma-1)));
    p2 = p1*(2*gamma*M1n^2-(gamma-1))/(gamma+1);
    rho2 = rho1*((gamma+1)*M1n^2)/((gamma-1)*M1n^2+2);
end