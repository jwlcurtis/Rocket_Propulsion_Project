function [T_a]=Air_Temprature(alt)
    if alt<32809
        T_a=-0.0036*alt+518;
    else
        T_a=399;
    end
end