function [mp]=Prop_mass(N,r1,r0,L0,wi,rho_p)
    mp=N*pi*(r0^2-(r1+wi)^2)*rho_p*(L0-2*wi);
end