% close all
clear all
clc
format long
global P
global a 
global colour_r1 colour_r2

tic

% numerical config
tfinal = 300*100; % final time in days
age_max = 50*365; % max ages in days
dt = 30; % time/age step size in days
da = dt;
t = (0:dt:tfinal)'; nt = length(t);
a = (0:da:age_max)'; na = length(a);

% model parameters - rates are in 1/day
baseline_Malaria_parameters;

% allocation
% SH, EH, etc.. = cell averages
SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt); RH = NaN(na,nt);
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt); 
Cs = NaN(na,nt);
rho_ave = NaN(1,nt);
NM = P.gM/P.muM;
NH = 1;

% initial condition 
[SH(:,1),EH(:,1),DH(:,1),AH(:,1),RH(:,1),SM(1,1),EM(1,1),IM(1,1)] = Malaria_IC(NH,NM); 
Cs(:,1) = Immunity_IC; % initial immunity and related probability
rho_ave(1,1) = mean(P.rho);

%% time evolution
for n = 1:nt-1
    NHa = SH(:,n)+EH(:,n)+DH(:,n)+AH(:,n)+RH(:,n); % total human at age a, t = n
    NH = trapz(NHa)*da; % total human population at t=n;
    
    NM = SM(1,n)+EM(1,n)+IM(1,n);
    [bH,~] = biting_rate(NH,NM); 
    lamH = FOI_H(bH,IM(1,n),NM);  % force of infection at t=n
    
    % human birth terms
    SH(1,n+1) = trapz(P.gH.*NHa)*da;
    EH(1,n+1) = 0;
    DH(1,n+1) = 0;
    AH(1,n+1) = 0;
    RH(1,n+1) = 0;
    % maternal immunity at age = 0
    Cs(1,n+1) = trapz(P.gH.*Cs(:,n).*NHa)/NH*da;
    
    for k = 1:na-1
        SH(k+1,n+1) = (SH(k,n)+dt*P.w(k+1)*RH(k,n)) / (1+(lamH+P.v(k+1)+P.muH(k+1))*dt);
        EH(k+1,n+1) = (EH(k,n)+dt*lamH*SH(k+1,n+1)) / (1+(P.h+P.muH(k+1))*dt);
        % solve DH and AH together
        num = (1+dt*(P.rA+ P.muH(k+1)))*(DH(k,n)+dt*P.h*P.rho(k+1)*EH(k+1,n+1))+...
            dt*(AH(k,n)+DH(k,n)+dt*P.h*EH(k+1,n+1))*P.psi(k+1)*lamH;
        den = (1+dt*(P.rA+ P.muH(k+1)))*(1+dt*(P.rD+P.muD(k+1)+P.muH(k+1)))+...
            dt*(1+dt*(P.muD(k+1)+P.muH(k+1)+P.rD*P.phi(k+1)))*P.psi(k+1)*lamH;
        DH(k+1,n+1) = num/den;
        AH(k+1,n+1) = (AH(k,n)+dt*((1-P.rho(k+1))*P.h*EH(k+1,n+1)+(1-P.phi(k+1))*P.rD*DH(k+1,n+1))) / (1+dt*(P.psi(k+1)*lamH+P.rA+P.muH(k+1)));
        RH(k+1,n+1) = ((1-P.w(k+1)*dt)*RH(k,n)+dt*(P.phi(k+1)*P.rD*DH(k+1,n+1)+P.rA*AH(k+1,n+1)+P.v(k+1)*SH(k+1,n+1)))/(1+dt*P.muH(k+1));
        % constraint on timestep for positivity: (1-P.w(k+1)*dt) > 0               
        
%         % update immunity use Qn
%         Qn = P.c1*lamH*SH(k,n)/NH + P.c2*lamH*AH(k,n)/NH;
%         Vn = P.c3*Cs(1,n+1).*exp(-a(k+1)/P.dm);
%         Cs(k+1,n+1) = (Cs(k,n)+dt*(Qn+(1/P.ds-1)*Vn))/(1+1/P.ds*dt); 
    end
        
    % adjust mosquito infection accordingly - use tn level!
    [SM(1,n+1),EM(1,n+1),IM(1,n+1)] = mosquito_ODE(DH(:,n),AH(:,n),NH,NM);
    
    % update immunity use Qn+1
    NHap1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1)+RH(:,n+1); % total human at age a, t = n
    NHp1 = trapz(NHap1)*da; % total human population at t=n;
    NMp1 = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    [bHp1,~] = biting_rate(NHp1,NMp1);
    lamHp1 = FOI_H(bHp1,IM(1,n+1),NMp1);   
    for k = 1:na-1       
        Qnp1 = P.c1*lamHp1*SH(k+1,n+1)/NHp1 + P.c2*lamHp1*AH(k+1,n+1)/NHp1;       
        Vn = P.c3*Cs(1,n+1).*exp(-a(k+1)/P.dm);
        Cs(k+1,n+1) = (Cs(k,n)+dt*(Qnp1+(1/P.ds-1)*Vn))/(1+1/P.ds*dt); % use Qn+1
    end
    
    % update progression probability based on immunity Cs
    P.phi = sigmoid_prob(Cs(:,n+1), 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Cs(:,n+1), 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Cs(:,n+1), 'psi'); % prob. AH -> DH
    rho_ave(1,n+1) = mean(P.rho);
end

%% Plotting
a_year = a/365;
age_group = NaN(age_max/365,2);
for i=1:age_max/365 % index for age i = age_group(i,1),..., age_group(i,2)
    temp = find((a_year>=i-1)&(a_year<i));
    age_group(i,1) = min(temp); 
    age_group(i,2) = max(temp);
end
toc;

% figure_setups;
% plot(t,trapz(SH,1)*da,'b-*'); hold on;
% plot(t,trapz(EH,1)*da,'-.','Color',colour_r1);
% plot(t,trapz(AH,1)*da,'-o','Color',colour_r2);
% plot(t,trapz(DH,1)*da,'r-');
% plot(t,trapz(RH,1)*da,'g-');
% plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1)+trapz(RH,1))*da)
% legend('SH-age','EH-age','AH-age', 'DH-age','RH-age','$N_H$');
% title('human population size by stages')
% axis_years(gca,tfinal); % change to x-axis to years if needed
% grid on
% % axis([0 tfinal 0 1.1])


figure_setups; 
Cs_year = NaN(age_max/365,nt); % cell average by ages
for n = nt
    for i = 1:age_max/365
        Cs_year(i,n) = mean(Cs(age_group(i,1):age_group(i,2),n),1);
    end
    plot(Cs_year(:,n));
    title(['~~~Immun dist at t = ', num2str(t(n)/365,3), 'yrs'])
%     if n==1; pause; else; pause;end
end
xlabel('age (years)')
% ylim([7.4 8.8]*10^-5)
grid on

figure_setups;
plot(t,trapz(Cs,1)*da);
% axis([0 tfinal 0 10])
axis_years(gca,tfinal)
title('~~~Total Immun in time');
grid on

figure_setups;
plot(a,P.rho)
axis_years(gca,tfinal)
title('$\rho$(age) at tfinal')
grid on

figure_setups;
plot(t,rho_ave)
axis_years(gca,tfinal)
title('average $\rho$ in time')
grid on


% figure_setups;
% plot(t,SM,'b-'); hold on;
% plot(t,EM,'-','Color',colour_r1);
% plot(t,IM,'r-');
% plot(t,SM+EM+IM)
% legend('SM','EM','IM','$N_M$');
% title('mosquito population size by stages')
% axis_years(gca,tfinal); % change to x-axis to years if needed
% grid on
% % axis([0 tfinal 0 5])

% %% final age distributions
% figure_setups;
% plot(a/365,SH(:,end),'b-'); hold on;
% plot(a/365,EH(:,end),'-','Color',colour_r1);
% plot(a/365,AH(:,end),'-','Color',colour_r2);
% plot(a/365,DH(:,end),'r-');
% plot(a/365,RH(:,end),'g-');
% legend('SH-age','EH-age','AH-age', 'DH-age','RH-age');
% title('Age distributions by class (final time)')
% grid on


