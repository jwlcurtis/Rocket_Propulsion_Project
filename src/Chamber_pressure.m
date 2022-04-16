function [P_c]=Chamber_pressure(a0,sigmap,Tb,Tb0,rho_p,c_star,Ab,g,At,n)
    P_c =(a0*exp(sigmap*(Tb-Tb0))*rho_p*c_star*Ab/(g*At))^(1/(1-n));
end