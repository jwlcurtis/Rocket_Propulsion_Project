clc
clear
format long
%% MAE 440 PR08D
% Curtis and Braylyan

%% User Changeable Variables
N = 4;
L0 = 3.00;      %in
r1 = 1.71;      %in
r0 = 2.375;     %in
At(1) = 1.5;   %in;
AeAt(1) = 1.5;
Tb = 30;      % F
m_ballast = 0.18;  %lbm
afterbo = true;
t_end = 40; % arbitrary ending time
%% Non-User Variables
g = 32.2;       % [ft/s^2]
a0 = 0.030;     % [in/s *psia^-n]
n = 0.35;
sigmap = 0.001;  % [/deg F]
gamma = 1.25;
c_star = 5210;  % [ft/s]
rho_p = 0.065;   % [lbm/in^3]
Tb0 = 70;        % [F]
m_payload = 40; % [lbm]
dw = 0.01;      % [in]
alt(1) = 0;     % [ft]
gamma_air = 1.4;% [unitless]
R_u = 1545.3;   % [lbf ft/ lb-mol R]
mw_air = 28.97;  % [lbm/lbmol]
Rocket_d=6.19;  % [in]
%% Burn Calculations 
Ae = AeAt(1)*At(1);   % Area of nozzle exit
d_throat(1) = 2*sqrt(At(1)/pi);
case_length=Case_Length(N,L0);
case_mass=Case_Mass(case_length);
A_Rocket=Rocket_area(Rocket_d);
w=Web(r0,r1,dw,L0);

% internal ballistics
for i=1:(length(w)-1)
    if i==1
        t(i) = 0;   % time elapsed
    end
    
    % Internal Balistics
    A_burn(i)=Burn_area(N,r1,r0,L0,w(i));
    Mass_p(i)=Prop_mass(N,r1,r0,L0,w(i),rho_p);
    P_chamber(i)=Chamber_pressure(a0,sigmap,Tb,Tb0,rho_p,c_star,A_burn(i),g,At(i),n);
    r_Burn(i)=Burn_rate(a0,sigmap,Tb,Tb0,P_chamber(i),n);
    t(i+1)=t(i)+((w(i+1)-w(i))/r_Burn(i));
    d_throat(i+1)=d_throat(i)+0.000087*(t(i+1)-t(i))*P_chamber(i);
    At(i+1)=(d_throat(i+1)/2)^2*pi;
    AeAt(i+1)=Ae/At(i+1);
end

A_burn(length(w))=Burn_area(N,r1,r0,L0,w(length(w)));
Mass_p((length(w)))=Prop_mass(N,r1,r0,L0,w((length(w))),rho_p);
P_chamber((length(w)))=Chamber_pressure(a0,sigmap,Tb,Tb0,rho_p,c_star,A_burn((length(w))),g,At((length(w))),n);
r_Burn((length(w)))=Burn_rate(a0,sigmap,Tb,Tb0,P_chamber((length(w))),n);

% cf and vechical dynamics
for i=1:(length(w)-1)
    if i==1
        t(i) = 0;   % time elapsed
        v_rocket(i) = 0;
        alt(i) =0;
    end
    % atmosphere calcs
    Temp_air(i)=Air_Temprature(alt(i));
    Pressure_air(i)=Air_Pressure(alt(i));
    Density_air(i)=Air_Density(alt(i));
    Pa_Pc(i)=Pressure_air(i)/P_chamber(i);
    
    %thrust
    cf(i)=Cf_solver(AeAt(i),Pa_Pc(i),2,gamma);
    F(i)=At(i)*P_chamber(i)*cf(i);
    
    % CD calc
    a_air(i)=Speed_of_sound_air(gamma_air,R_u,mw_air,g,Temp_air(i));
    Rocket_M(i)=Rocket_Mach(v_rocket(i),a_air(i));
    CD(i)=Drag_Coefficent(Rocket_M(i));
    % vechical dynamics
    m_Rocket(i) = Mass_p(i)+m_ballast+case_mass+m_payload;
    D(i) = 1/2*Density_air(i)*CD(i)*v_rocket(i)^2*A_Rocket; % Drag
    FM(i) = F(i)/m_Rocket(i) * 32.2;
    DM(i) = D(i)/m_Rocket(i);
    accl(i) = FM(i) - DM(i) - g ;                % accl in ft/s^2
    accl_g(i) = accl(i)/g;% accl in g's
    v_rocket(i+1) = accl(i)*(t(i+1)-t(i))+v_rocket(i);
    alt(i+1) = alt(i)+(v_rocket(i)+v_rocket(i+1))/2*(t(i+1)-t(i));
end
    Temp_air(length(w))=Air_Temprature(alt(length(w)));
    Pressure_air((length(w)))=Air_Pressure(alt(length(w)));
    Density_air((length(w)))=Air_Density(alt(length(w)));
    Pa_Pc((length(w)))=Pressure_air((length(w)))/P_chamber((length(w)));
    cf((length(w)))=Cf_solver(AeAt((length(w))),Pa_Pc((length(w))),3,gamma);
    F((length(w)))=At((length(w)))*P_chamber((length(w)))*cf((length(w)));
    
    % CD calc
    a_air((length(w)))=Speed_of_sound_air(gamma_air,R_u,mw_air,g,Temp_air((length(w))));
    Rocket_M((length(w)))=Rocket_Mach(v_rocket((length(w))),a_air((length(w))));
    CD((length(w)))=Drag_Coefficent(Rocket_M((length(w))));
    % vechical dynamics
    m_Rocket((length(w))) = Mass_p((length(w)))+m_ballast+case_mass+m_payload;
    D((length(w))) = 1/2*Density_air((length(w)))*CD((length(w)))*v_rocket((length(w)))^2*A_Rocket; % Drag
    FM((length(w))) = F((length(w)))/m_Rocket((length(w))) * 32.2;
    DM((length(w))) = D((length(w)))/m_Rocket((length(w)));
    accl((length(w))) = FM((length(w))) - DM((length(w))) - g ;                % accl in ft/s^2
    accl_g((length(w))) = accl((length(w)))/g;% accl in g's
    t_burnout=t(end);
 
%% impulse
for i=1:length(F)
    if i==1
        I(i) = 0;
    else
        I(i) = (F(i)+F(i-1))/2*(t(i)-t(i-1));
    end
end

I_total=sum(I)
ISP=I_total/Mass_p(1);
%% Coast Phase
if afterbo == true
    dt = 0.1;
    
    m_Rocket_final=m_Rocket(end);
    while t<(t_end-dt)
        i = i+1;
        t(i) = t(end)+dt;
        Mass_p(i) = 0;
        P_chamber(i) = 0;
        v_rocket(i) = accl(i-1)*(t(i)-t(i-1))+v_rocket(i-1);
        alt(i) = alt(i-1)+(v_rocket(i-1)+v_rocket(i))/2*(t(i)-t(i-1));
        
        Temp_air(i)=Air_Temprature(alt(i));
        Pressure_air(i)=Air_Pressure(alt(i));
        Density_air(i)=Air_Density(alt(i));
        
        
        % Speed of sound
        a_air(i)=Speed_of_sound_air(gamma_air,R_u,mw_air,g,Temp_air(i));
        
        % Mach
        Rocket_M(i)=Rocket_Mach(v_rocket(i),a_air(i));
        if Rocket_M(i)<0.6 
            CD(i)=0.15;
        end
        if Rocket_M(i)>=0.6 && Rocket_M(i)<1.2
            CD(i)=-0.12+0.45*Rocket_M(i);
        end
        if Rocket_M(i)>=1.2 && Rocket_M(i)<=1.8
            CD(i)=0.76-0.283*Rocket_M(i);
        end
        if Rocket_M(i)>1.8 && Rocket_M(i)<=4.0
            CD(i)=0.311-0.034*Rocket_M(i);
        end
        if Rocket_M(i)>4.0
            CD(i)=0.175;
        end
        
        
        m_Rocket(i) = m_Rocket_final;
        D(i) = 1/2*Density_air(i)*CD(i)*v_rocket(i)^2*A_Rocket; % Drag
        DM(i) = D(i)/m_Rocket(i);
        accl(i) = -DM(i) - g ;                % accl in ft/s^2
        accl_g(i) = accl(i)/g;                      % accl in g's
        F(i) = 0;
    end
end

%% Requested Values
% Table 1-1 values
T1_1 = table(a0, n, sigmap, gamma, c_star, rho_p, mw_air, gamma_air,...
    A_Rocket, R_u)

% Table 2 values
T1_2 = table(r1, r0, L0, N, AeAt(1),At(1),m_ballast)

% Table 1-3 values
At_bo = At(end);
ApAt = pi*r1^2/At(1);
F_max=max(F);
Pc_max=max(P_chamber);

mp_total=Mass_p(1);
Ap0_At0=pi*r1^2/At(1);
T1_3 = table(mp_total, At_bo, t_burnout, ISP,Pc_max,F_max,Ap0_At0)

% Table 1-4 values
m0=m_Rocket(1);
h_bo=alt(length(w));
v_p=v_rocket(length(w));
a_p=accl(length(w));
h_max=max(alt);
v_max=max(v_rocket);
a_max=max(accl);
g_max=a_max/g;
T1_3 = table(case_mass, m0 , h_bo, v_p,a_p,h_max,v_max,a_max,g_max)
%% Graphs
figure (1)
% ti1=t(length(t)):0.1:5;
% for i=1:length(ti1)
%     mpi1(i)=0;
% end
% ti2=cat(2,t,ti1);
% mpi=cat(2,Mass_p,mpi1);
plot(t,Mass_p,'LineWidth',1.5)
xlabel("Time [s]")
ylabel("Mass of Propellent [lbs]")

figure (2)
% for i=1:length(ti1)
%     pci1(i)=0;
% end
% pci=cat(2,P_chamber,pci1);
plot(t,P_chamber,'LineWidth',1.5)
xlabel("Time [s]")
ylabel("Chamber Pressure [psi]")

figure (3)

% for i=1:length(ti1)
%     F1(i)=0;
% end
% F=cat(2,F,F1);
plot(t,F,'LineWidth',1.5)
xlabel("Time [s]")
ylabel("Thrust [lbf]")

figure (4)
% ttotal=cat(2,t,tcoast);
% Atotal=cat(2,accl,accl_coast);
plot(t,accl,'LineWidth',1.5)
xlabel("Time [s]")
ylabel("Accelration [ft/s^2]")

figure (5)
% Vtotal=cat(2,v_rocket,v_rocket_coast);
plot(t,v_rocket,'LineWidth',1.5)
xlabel("Time [s]")
ylabel("Velocity [ft/s]")

figure (6)
% Dtotal=cat(2,alt,alt_coast);
plot(t,alt,'LineWidth',1.5)
xlabel("Time [s]")
ylabel("Altitude [ft]")

Cf_graph(AeAt,cf,gamma,7)
%figure(8)
%plot(AeAt,t)
%% Functions

% % Table 1-4 values
% max_alt = max(alt);
% max_v = max(v_rocket);
% accl_max = max(accl);
% accl_max_g = max(accl_g);
% T1_4 = table(prop_case_mass, MTOW,h_bo,v_bo,a_bo,max_alt,max_v,...
%     accl_max,accl_max_g)

% function [Ab]=Burn_area(N,r1,r0,L0,wi)
%     Ab=N*2*pi*((r1+wi)*(L0-2*wi)+(r0^2-(r1+wi)^2));
% end
% 
% function [mp]=Prop_mass(N,r1,r0,L0,wi,rho_p)
%     mp=N*pi*(r0^2-(r1+wi)^2)*rho_p*(L0-2*wi);
% end
% function [w]=Web(r0,r1,dw,L0)
%     if (L0/2) < (r0-r1)
%     w_max = L0/2;
%     else
%     w_max = r0-r1;
%     end
%     w=0:dw:w_max;
%     if w(length(w))~=w_max
%     w(length(w)+1)=w_max;
%     end
%     
% end
% 
% function [A_rocket]=Rocket_area(Rocket_d)
%     A_rocket=(pi*(Rocket_d/2)^2)/144; %[ft^2]
% end
% 
% function [P_a]= Air_Pressure(alt)
%     if alt<83000
%         P_a=(-4.272981*10^(-14))*alt^3+(0.000000008060081*alt^2)-0.0005482655*alt + 14.69241;
%     else
%         P_a=0;
%     end
% end
% 
% function [T_a]=Air_Temprature(alt)
%     if alt<32809
%         T_a=-0.0036*alt+518;
%     else
%         T_a=399;
%     end
% end
% 
% function [rho_a]=Air_Density(alt)
%     if alt<82000
%         rho_a=0.00000000001255*alt^2-(0.0000019453*alt)+0.07579;
%     else
%         rho_a=0;
%     end
% end
% 
% function [cd1]=Drag_Coefficent(mach)
%     if mach<0.6 && mach>=0
%         cd1=0.15;
%     end
%     if mach>=0.6 && mach<1.2
%         cd1=-0.12+0.45*mach;
%     end
%     if mach>=1.2 && mach<=1.8
%         cd1=0.76-0.283*mach;
%     end
%     if mach>1.8 && mach<=4.0
%         cd1=0.311-0.034*mach;
%     end
%     if mach>4.0
%         cd1=0.175;
%     end
% end
% 
% function [a_air]=Speed_of_sound_air(gamma_air,R_u,mw_air,g,T_a)
%     a_air=sqrt(gamma_air*(R_u/mw_air)*T_a*g);
% end
% 
% function [Mach_R]=Rocket_Mach(Rocket_V,a_air)
%     Mach_R=Rocket_V/a_air;
% end
% 
% function [Mach_E]=Engine_Mach(e,gamma,m0)
%     mach=@(m) (1/m)*((2+(gamma-1)*m^2)/(gamma+1))^((gamma+1)/(2*(gamma-1)))-(e);
%     Mach_E=fzero(mach,m0);
% end
% 
% function [L_case]=Case_Length(N,L0)
%     L_case=0.125*(N)+(N*L0);
% end
% 
% function [m_case]=Case_Mass(L_case)
%     m_case=0.25*L_case;
% end
% 
% function []=Cf_graph(gamma,fig)
%     e=1:0.01:100;
%     cf2=zeros(1,length(e));
%     cf4=zeros(1,length(e));
%     cf10=zeros(1,length(e));
%     cf25=zeros(1,length(e));
%     cf68=zeros(1,length(e));
%     cf150=zeros(1,length(e));
%     cf500=zeros(1,length(e));
%     cf1000=zeros(1,length(e));
%     cfinf=zeros(1,length(e));
%     cfmin=zeros(1,length(e));
%     cfopt=zeros(1,length(e));
%     for i=1:length(e)
%         if i==1
%             m0=1;
%         else
%         m0=2; % inital guess
%         end
%         cf2(i)=Cf_solver(e(i),1/2,m0,gamma);
%         cf4(i)=Cf_solver(e(i),1/4,m0,gamma);
%         cf10(i)=Cf_solver(e(i),1/10,m0,gamma);
%         cf25(i)=Cf_solver(e(i),1/25,m0,gamma);
%         cf68(i)=Cf_solver(e(i),1/68,m0,gamma);
%         cf150(i)=Cf_solver(e(i),1/150,m0,gamma);
%         cf500(i)=Cf_solver(e(i),1/500,m0,gamma);
%         cf1000(i)=Cf_solver(e(i),1/1000,m0,gamma);
%         cfinf(i)=Cf_solver(e(i),0,m0,gamma);
%         cfmin(i)=-0.0445*log(e(i))^2+0.5324*(log(e(i)))+0.1843;
%         cfopt(i)=Cf_opt_solver(e(i),0,m0,gamma);
%     end
%     figure(fig)
%     semilogx(e,cf2)
%     hold on
%     semilogx(e,cf4)
%     semilogx(e,cf10)
%     semilogx(e,cf25)
%     semilogx(e,cf68)
%     semilogx(e,cf150)
%     semilogx(e,cf500)
%     semilogx(e,cf1000)
%     semilogx(e,cfinf)
%     semilogx(e,cfmin,"--")
%     semilogx(e,cfopt,"--")
%     hold off
%     ylabel("c_f")
%     xlabel("A_e/A_t")
%     ylim([0.6 2])
%     xlim([1 100])
%     legend({"P_c/P_a=2","P_c/P_a=4","P_c/P_a=10","P_c/P_a=25","P_c/P_a=68","P_c/P_a=150","P_c/P_a=500","P_c/P_a=1000","P_c/P_a=\infty","Line of Seperation","Optimal c_f"},'Location','northwest','NumColumns',3)
% end
% 
% function [cf]=Cf_solver(e,Pa_Pc,m0,gamma)
%     Mach_E=Engine_Mach(e,gamma,m0);
%     Pc_Pe=(1+((gamma-1)/2)*(Mach_E^2))^(gamma/(gamma-1));
%     papc=Pa_Pc;
%     pepc=1/Pc_Pe;
%     exp1=(gamma-1)/gamma;
%     cf=sqrt((((2*gamma^2)/(gamma-1))*((2/(gamma+1))^((gamma+1)/(gamma-1)))*((1-((pepc)^(exp1))))))+(e*((pepc-papc)));
% end
% 
% function [cf]=Cf_opt_solver(e,Pa_Pc,m0,gamma)
%     Mach_E=Engine_Mach(e,gamma,m0);
%     Pc_Pe=(1+((gamma-1)/2)*(Mach_E^2))^(gamma/(gamma-1));
%     papc=Pa_Pc;
%     pepc=1/Pc_Pe;
%     exp1=(gamma-1)/gamma;
%     cf=sqrt((((2*gamma^2)/(gamma-1))*((2/(gamma+1))^((gamma+1)/(gamma-1)))*((1-((pepc)^(exp1))))));
% end
% 
% function [P_c]=Chamber_pressure(a0,sigmap,Tb,Tb0,rho_p,c_star,Ab,g,At,n)
%     P_c =(a0*exp(sigmap*(Tb-Tb0))*rho_p*c_star*Ab/(g*At))^(1/(1-n));
% end
% 
% function [r_b]=Burn_rate(a0,sigmap,Tb,Tb0,P_c,n)
%     r_b=a0*(exp(sigmap*(Tb-Tb0)))*(P_c^n);
% end
