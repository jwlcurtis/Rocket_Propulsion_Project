function [r_b]=Burn_rate(a0,sigmap,Tb,Tb0,P_c,n)
    r_b=a0*exp(sigmap*(Tb-Tb0))*P_c^n;
end