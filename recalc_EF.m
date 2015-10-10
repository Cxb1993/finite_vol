function [E,F,rho,u,v,p] = recalc_EF(U,E,F)
        %Get bar values for flux terms E and F from U bar
            rho = U(:,:,1);
            u = U(:,:,2)./rho;
            v = U(:,:,3)./rho;
            e = U(:,:,4);
            p = pressure(rho,e,u,v); 
        %E flux
            E(:,:,1) = rho.*u;
            E(:,:,2) = rho.*u.^2+p;
            E(:,:,3) = rho.*u.*v;
            E(:,:,4) = (e+p).*u;
        %F flux
            F(:,:,1) = rho.*v;
            F(:,:,2) = E(:,:,3);
            F(:,:,3) = rho.*v.^2+p;
            F(:,:,4) = (e+p).*v;
%--------------------------------------------------------------------------
%Auxillarily Functions
%--------------------------------------------------------------------------
    function [p] = pressure(rho, e, u,v)
    gamma = 1.4;
    p = (gamma-1)*(e-rho.*(u.^2+v.^2)/2);
    end
end