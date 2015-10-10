function [beta] = tbM(theta,M)
%Initial Beta
    beta = 10*pi/180;
%use newtons method to solve for beta
    db = 10^-8;
    converge = 1;
    while converge > 0.00001
        tant = tantheta(M,beta);
        tant_prime = (tantheta(M,beta+db)-tantheta(M,beta-db))/2/db;
        beta2 = beta-(tant-tan(theta*pi/180))/tant_prime;
        converge = abs(beta2-beta);
        beta = beta2;
    end
    beta = beta*180/pi;
end