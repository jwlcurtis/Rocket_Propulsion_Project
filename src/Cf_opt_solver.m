function [cf]=Cf_opt_solver(e,Pa_Pc,m0,gamma)
    Mach_E=Engine_Mach(e,gamma,m0);
    Pc_Pe=(1+((gamma-1)/2)*(Mach_E^2))^(gamma/(gamma-1));
    papc=Pa_Pc;
    pepc=1/Pc_Pe;
    exp1=(gamma-1)/gamma;
    cf=sqrt((((2*gamma^2)/(gamma-1))*((2/(gamma+1))^((gamma+1)/(gamma-1)))*((1-((pepc)^(exp1))))));
end
function [P_c]=Chamber_pressure(a0,sigmap,Tb,Tb0,rho_p,c_star,Ab,g,At,n)
    P_c =(a0*exp(sigmap*(Tb-Tb0))*rho_p*c_star*Ab/(g*At))^(1/(1-n));
end