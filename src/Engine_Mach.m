function [Mach_E]=Engine_Mach(e,gamma,m0)
    mach=@(m) (1/m)*((2+(gamma-1)*m^2)/(gamma+1))^((gamma+1)/(2*(gamma-1)))-(e);
    Mach_E=fzero(mach,m0);
end