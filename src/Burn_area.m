function [Ab]=Burn_area(N,r1,r0,L0,wi)
    Ab = N*(2*pi*(r1+wi)*(L0-2*wi))+...
            N*(2*pi*(r0^2-(r1+wi)^2));
end
