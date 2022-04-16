function [P_a]= Air_Pressure(alt)
    if alt<83000
        P_a=(-4.272981*10^(-14))*alt^3+(0.000000008060081*alt^2)-0.0005482655*alt + 14.69241;
    else
        P_a=0;
    end
end