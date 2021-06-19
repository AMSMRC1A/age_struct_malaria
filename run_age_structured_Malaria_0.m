%% run codes with big matrices
% close all
clear all
clc
format long
global P
global colour_r1 colour_r2

tic

% numerical config
tfinal = 50*365; % final time in days
age_max = 50*365; % max ages in days
P.age_max = age_max;
dt = 365; % time/age step size in days
da = dt;
t = (0:dt:tfinal)'; nt = length(t);
a = (0:da:age_max)'; na = length(a);

P.a = a;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;

% model parameters - rates are in 1/day
baseline_Malaria_parameters;

% allocation
% SH, EH, etc.. = cell averages
SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt); RH = NaN(na,nt);
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt); 
Cs = NaN(na,nt); Cm = NaN(na,nt); Cac = NaN(na,nt); 
rho_ave = NaN(1,nt);
NM = P.gM/P.muM;
NH = 1;

% initial condition 
[SH(:,1),EH(:,1),DH(:,1),AH(:,1),RH(:,1),SM(1,1),EM(1,1),IM(1,1)] = Malaria_IC(NH,NM); 
[Cm(:,1),Cac(:,1),Cs(:,1)] = Immunity_IC; % initial immunity and related probability
rho_ave(1,1) = mean(P.rho);

%% time evolution
for n = 1:nt-1
    if mod(n,(nt-1)/10)==0
        display(['progress = ',num2str(n/(nt-1)*100),'%'])
    end
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
    
    SH(2:end,n+1) = (SH(1:end-1,n)+P.dt*P.w(2:end).*RH(1:end-1,n)) ./ (1+(lamH+P.v(2:end)+P.muH(2:end))*P.dt);
    EH(2:end,n+1) = (EH(1:end-1,n)+P.dt*lamH*SH(2:end,n+1)) ./ (1+(P.h+P.muH(2:end))*P.dt);
    % solve DH and AH together
        num = (1+P.dt*(P.rA+ P.muH(2:end))).*(DH(1:end-1,n)+P.dt*P.h*P.rho(2:end).*EH(2:end,n+1))+...
            P.dt*(AH(1:end-1,n)+DH(1:end-1,n)+P.dt*P.h*EH(2:end,n+1)).*P.psi(2:end)*lamH;
        den = (1+P.dt*(P.rA+ P.muH(2:end))).*(1+P.dt*(P.rD+P.muD(2:end)+P.muH(2:end)))+...
            P.dt*(1+P.dt*(P.muD(2:end)+P.muH(2:end)+P.rD*P.phi(2:end))).*P.psi(2:end)*lamH;
    DH(2:end,n+1) = num./den;
	AH(2:end,n+1) = (AH(1:end-1,n)+P.dt*((1-P.rho(2:end))*P.h.*EH(2:end,n+1)+(1-P.phi(2:end))*P.rD.*DH(2:end,n+1))) ./ (1+P.dt*(P.psi(2:end)*lamH+P.rA+P.muH(2:end)));
	RH(2:end,n+1) = ((1-P.w(2:end)*P.dt).*RH(1:end-1,n)+P.dt*(P.phi(2:end)*P.rD.*DH(2:end,n+1)+P.rA*AH(2:end,n+1)+P.v(2:end).*SH(2:end,n+1))) ./ (1+P.dt*P.muH(2:end));
        
    % adjust mosquito infection accordingly - use tn level!
    [SM(1,n+1),EM(1,n+1),IM(1,n+1)] = mosquito_ODE(DH(:,n),AH(:,n),NH,NM);
    
    % immunity at age = 0
    Cm(1,n+1) = P.m*trapz(P.gH.*Cs(:,n).*NHa)/NH*da; % why is Cs in the integral here?
    Cac(1,n+1) = 0;
    % maternal immunity
    n0 = min(n,na-1);
    Cm(2:n0+1,n+1) = (Cm(1,1:n0))'.*exp(-a(2:n0+1)/P.dm); % k=1:n0  
    Cm(n0+2:end,n+1) = Cm(2:end-n0,1).*exp(-t(n+1)/P.dm);  % k=n0+1:end-1   
    % acquired immunity
    NHap1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1)+RH(:,n+1); % total human at age a, t = n
    NHp1 = trapz(NHap1)*da; % total human population at t=n;
    NMp1 = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    [bHp1,~] = biting_rate(NHp1,NMp1);
    lamHp1 = FOI_H(bHp1,IM(1,n+1),NMp1);   
    Qnp1 = lamHp1*(P.w1*SH(2:end,n+1) + P.w2*EH(2:end,n+1) + P.w3*AH(2:end,n+1) ...
        + P.w4*DH(2:end,n+1) + P.w5*RH(2:end,n+1))/NHp1;
    Cac(2:end,n+1) = (Cac(1:end-1,n)+P.dt*Qnp1)./(1+1/P.ds*P.dt); % use Qn+1
    Cs(:,n+1) = P.c1*Cac(:,n+1)+P.c2*Cm(:,n+1); % total immunity from acquired and maternal sources
    
    % update progression probability based on immunity Cs
    P.phi = sigmoid_prob(Cs(:,n+1), 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Cs(:,n+1), 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Cs(:,n+1), 'psi'); % prob. AH -> DH
    rho_ave(1,n+1) = mean(P.rho);
end
%% Plotting
% a_year = a/365;
% age_group = NaN(age_max/365,2);
% for i=1:age_max/365 % index for age i = age_group(i,1),..., age_group(i,2)
%     temp = find((a_year>=i-1)&(a_year<i));
%     age_group(i,1) = min(temp); 
%     age_group(i,2) = max(temp);
% end
toc

% figure_setups;
% plot(t,trapz(SH,1)*da,'b-'); hold on;
% plot(t,trapz(EH,1)*da,'-','Color',colour_r1);
% plot(t,trapz(AH,1)*da,'-','Color',colour_r2);
% plot(t,trapz(DH,1)*da,'r-');
% plot(t,trapz(RH,1)*da,'g-');
% plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1)+trapz(RH,1))*da)
% legend('SH-age','EH-age','AH-age', 'DH-age','RH-age','$N_H$');
% title('human population size by stages')
% axis_years(gca,tfinal); % change to x-axis to years if needed
% grid on
% axis([0 tfinal 0 1.1])

% figure_setups; 
% Cs_year = NaN(age_max/365,nt); % cell average by yearly ages
% for n = nt
%     for i = 1:age_max/365
%         Cs_year(i,n) = mean(Cs(age_group(i,1):age_group(i,2),n),1);
%     end
%     plot(Cs_year(:,n));
%     title(['~~~Immun dist at t = ', num2str(t(n)/365,3), 'yrs'])
% %     if n==1; pause; else; pause;end
% end
% xlabel('age (years)')
% % ylim([7.4 8.8]*10^-5)
% grid on

%% Immunity related figures
% figure_setups; 
% subplot(2,2,1), plot(a,Cs(:,end));
% xlabel('age (years)')
% title(['~~~Final Immun dist, dt=', num2str(dt)])
% grid on
% 
% subplot(2,2,2), plot(t,trapz(Cs,1)*da);
% % axis_years(gca,tfinal)
% title(['~~~Pop. Immun vs time, dt=',num2str(dt)]);
% xlabel('time (days)')
% grid on
% 
% subplot(2,2,3), imagesc(t/365,a/365,Cs);
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('time');
% title('Total Immunity');

%% Impact of immunity on the sigmoids (rho, psi, phi)
figure_setups;
subplot(2,2,1), imagesc(t/365,a/365,sigmoid_prob(Cs, 'phi'));
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('$\phi$');

subplot(2,2,2), imagesc(t/365,a/365,sigmoid_prob(Cs, 'rho'));
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('$\rho$');

subplot(2,2,3), imagesc(t/365,a/365,sigmoid_prob(Cs, 'psi'));
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('$\psi$');

figure_setups;
plot(a,P.rho)
axis_years(gca,tfinal)
title('$\rho$(age) at tfinal')
grid on



%% Mosquito infection dynamics
% figure_setups;
% plot(t,SM,'b-'); hold on;
% plot(t,EM,'-','Color',colour_r1);
% plot(t,IM,'r-');
% plot(t,SM+EM+IM)
% legend('SM','EM','IM','$N_M$');
% title('mosquito population size by stages')
% axis_years(gca,tfinal); % change to x-axis to years if needed
% grid on
% axis([0 tfinal 0 5])

%% final age distributions
% figure_setups;
% plot(a/365,SH(:,end),'b-'); hold on;
% plot(a/365,EH(:,end),'-','Color',colour_r1);
% plot(a/365,AH(:,end),'-','Color',colour_r2);
% plot(a/365,DH(:,end),'r-');
% plot(a/365,RH(:,end),'g-');
% legend('SH-age','EH-age','AH-age', 'DH-age','RH-age');
% title('Age distributions by class (final time)')
% grid on

%% Compare initial and final age distribution
% figure_setups;
% Total_pop = SH+EH+AH+DH+RH;
% plot(a/365, Total_pop(:,end),'g-'); hold on;
% %plot(a/365,Total_pop(:,ceil(end/2)),'-','Color',colour_r1); 
% plot(a/365, P.n_tilde,'-','Color',colour_r2);
% plot(a/365, P.n_tilde*exp(P.p_hat*tfinal) ,'-.');
% legend('Final Age dist. (sim.)','Stable Age Dist. (initial)','Final Age dist. (theory)');
% title('Population Age Distributions')


