function [rho_a]=Air_Density(alt)
    if alt<82000
        rho_a=0.00000000001255*alt^2-(0.0000019453*alt)+0.07579;
    else
        rho_a=0;
    end
end