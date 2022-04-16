function [a_air]=Speed_of_sound_air(gamma_air,R_u,mw_air,g,T_a)
    a_air=sqrt(gamma_air*(R_u/mw_air)*T_a*g);
end