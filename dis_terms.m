function [dphi,dmhi,dphj,dmhj] = dis_terms(u,v,c,p,point,err)
    j = point(1); i = point(2);
        u_c_ph = 0.5*(u(j,i+1)+c(j,i+1)+u(j,i)+c(j,i));
        u_c_mh = 0.5*(u(j,i-1)+c(j,i-1)+u(j,i)+c(j,i));
        abs_p_term = abs(p(j,i+1)-2*p(j,i)+p(j,i-1));
        p_term = p(j,i+1)+2*p(j,i)+p(j,i-1);
    dphi = err*u_c_ph*(abs_p_term/p_term);
    dmhi = err*u_c_mh*(abs_p_term/p_term);
        u_c_ph = 0.5*(v(j+1,i)+c(j+1,i)+v(j,i)+c(j,i));
        u_c_mh = 0.5*(v(j-1,i)+c(j-1,i)+v(j,i)+c(j,i));
        abs_p_term = abs(p(j+1,i)-2*p(j,i)+p(j-1,i));
        p_term = p(j+1,i)+2*p(j,i)+p(j-1,i);
    dphj = err*u_c_ph*(abs_p_term/p_term);
    dmhj = err*u_c_mh*(abs_p_term/p_term);
end