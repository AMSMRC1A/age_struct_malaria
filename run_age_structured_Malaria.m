% close all
clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic

% numerical config
P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero
tfinal = 500*365; % final time in days
age_max = 60*365; % max ages in days
P.age_max = age_max;
dt = 25; % time/age step size in days
da = dt;
t = (0:dt:tfinal)'; 
nt = length(t);
a = (0:da:age_max)'; 
na = length(a);

P.a = a;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;

% model parameters - rates are in 1/day
baseline_Malaria_parameters;

% allocation
% SH, EH, etc.. = cell averages
SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt);
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt); 

% per-person immunity level or \tilde(Ctotal) \tilde(Cm) \tilde(Cac) on overleaf
Ctot = NaN(na,nt); Cm = NaN(na,nt); Cac = NaN(na,nt); 
NM = P.gM/P.muM;
NH = ones(1,length(t));

% initial condition 
[SH(:,1),EH(:,1),DH(:,1),AH(:,1),SM(1,1),EM(1,1),IM(1,1)] = Malaria_IC(NH(1),NM); 
[Cm(:,1),Cac(:,1),Ctot(:,1)] = Immunity_IC; % initial immunity and related probability

%% time evolution
for n = 1:nt-1
    if mod(n,(nt-1)/10)==0
        disp(['progress = ',num2str(n/(nt-1)*100),'%']);
    end
    PH = SH(:,n)+EH(:,n)+DH(:,n)+AH(:,n); % total human at age a, t = n
    NH(n) = trapz(PH)*da; % total human population at t=n;    
    NM = SM(1,n)+EM(1,n)+IM(1,n);
    [bH,~] = biting_rate(NH(n),NM); 
    lamH = FOI_H(bH,IM(1,n));  % force of infection at t=n
    
    % human birth terms
    SH(1,n+1) = trapz(P.gH.*PH)*da;
    EH(1,n+1) = 0;
    DH(1,n+1) = 0;
    AH(1,n+1) = 0;
    
    %disp(['n = ',num2str(n)])
    
    SH(2:end,n+1) = (SH(1:end-1,n)+P.dt*(P.phi(1:end-1)*P.rD.*DH(1:end-1,n)+P.rA*AH(1:end-1,n))) ./ (1+(lamH+P.muH(2:end))*P.dt);   
    EH(2:end,n+1) = (EH(1:end-1,n)+P.dt*lamH*SH(2:end,n+1)) ./ (1+(P.h+P.muH(2:end))*P.dt);
    AH(2:end,n+1) = ((1-P.dt*P.rA)*AH(1:end-1,n)+P.dt*((1-P.rho(1:end-1))*P.h.*EH(2:end,n+1)+(1-P.phi(1:end-1)).*P.rD.*DH(1:end-1,n)))...
        ./ (1+P.dt*(P.psi(1:end-1)*lamH+P.muH(2:end)));
	DH(2:end,n+1) = ((1-P.dt*P.rD)*DH(1:end-1,n)+P.dt*(P.rho(1:end-1)*P.h.*EH(2:end,n+1)+P.psi(1:end-1).*lamH.*AH(2:end,n+1)))...
        ./ (1+P.dt*(P.muH(2:end)+P.muD(2:end)));	
    
    % adjust mosquito infection accordingly - use tn level!
    [SM(1,n+1),EM(1,n+1),IM(1,n+1)] = mosquito_ODE(DH(:,n),AH(:,n),NH(n),NM);
    
    % immunity gained at age = 0
    Cm(1,n+1) = P.m*trapz(P.gH.*PH.*Cac(:,n))*da/SH(1,n+1);
    Cac(1,n+1) = 0;
    % maternal immunity
    n0 = min(n,na-1);
    Cm(2:n0+1,n+1) = (Cm(1,1:n0))'.*exp(-a(2:n0+1)/P.dm); % k=1:n0  
    Cm(n0+2:end,n+1) = Cm(2:end-n0,1).*exp(-t(n+1)/P.dm);  % k=n0+1:end-1   
    % acquired immunity - use Qn+1
    PHp1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1); % total human at age a, t = n+1
    NHp1 = trapz(PHp1)*da; % total human population at t=n;
    NMp1 = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    [bHp1,~] = biting_rate(NHp1,NMp1);
    lamHp1 = FOI_H(bHp1,IM(1,n+1));   
    % neglecting disease induced mortality in Cac
    Qnp1 = f(lamHp1).*(P.cS*SH(2:end,n+1) + P.cE*EH(2:end,n+1) + P.cA*AH(2:end,n+1) ...
        + P.cD*DH(2:end,n+1)) + P.cV*P.v(2:end).*SH(2:end,n+1) ;
    Cac(2:end,n+1) = (Cac(1:end-1,n)+P.dt*Qnp1)./(1 + P.dt*(1./P.dac + P.muD(2:end)));
    % Cm is now per person but Cac is pooled so need to multiply Cm by PH
    % to get total immunity contribution
    Ctot(:,n+1) = P.c1*Cac(:,n+1)+P.c2*Cm(:,n+1).*PHp1; % total immunity from acquired and maternal sources
    
    % update progression probability based on immunity Ctot
    P.phi = sigmoid_prob(Ctot(:,n+1), 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Ctot(:,n+1), 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Ctot(:,n+1), 'psi'); % prob. AH -> DH

end
PH = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end); % total human at age a, t = n
NH(end) = trapz(PH)*da;
%% Plotting
% a_year = a/365;
% age_group = NaN(age_max/365,2);
% for i=1:age_max/365 % index for age i = age_group(i,1),..., age_group(i,2)
%     temp = find((a_year>=i-1)&(a_year<i));
%     age_group(i,1) = min(temp); 
%     age_group(i,2) = max(temp);
% end
toc
%% 
figure_setups;
Nh = (trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da;
plot(t,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
plot(t,trapz(EH,1)*da,'--','Color',colour_mat3);
plot(t,trapz(AH,1)*da,'-.','Color',colour_mat2);
plot(t,trapz(DH,1)*da,'-','Color',colour_mat7);
plot(t,Nh,'k-')
legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$','Location','NorthWest');
title('Human population sizes');
axis_years(gca,tfinal); % change to x-axis to years if needed
grid on
axis([0 tfinal 0 max(Nh)+0.1]);
%%
figure_setups;
plot(t,trapz(SH,1)*da./Nh,'-','Color',colour_mat1); hold on;
plot(t,trapz(EH,1)*da./Nh,'--','Color',colour_mat3);
plot(t,trapz(AH,1)*da./Nh,'-.','Color',colour_mat2);
plot(t,trapz(DH,1)*da./Nh,'-','Color',colour_mat7);
plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da./Nh,'k-')
legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$');
title('human population proportions');
axis_years(gca,tfinal); % change to x-axis to years if needed
xlabel('time');
grid on
axis([0 tfinal 0 1.1]);

% figure_setups; 
% Ctot_year = NaN(age_max/365,nt); % cell average by yearly ages
% for n = nt
%     for i = 1:age_max/365
%         Ctot_year(i,n) = mean(Ctot(age_group(i,1):age_group(i,2),n),1);
%     end
%     plot(Ctot_year(:,n));
%     title(['~~~Immun dist at t = ', num2str(t(n)/365,3), 'yrs'])
% %     if n==1; pause; else; pause;end
% end
% xlabel('age (years)')
% % ylim([7.4 8.8]*10^-5)
% grid on

%% Immunity related figures
figure_setups; 
PH = SH+EH+DH+AH;
subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/4))./PH(:,floor(nt/4)));
hold on;
subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/2))./PH(:,floor(nt/2)));
subplot(2,2,1), plot(a/365,Ctot(:,floor(3*nt/4))./PH(:,floor(3*nt/4)));
subplot(2,2,1), plot(a/365,Ctot(:,end)./PH(:,end));
xlabel('age')
legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
    ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)],'Location','NorthWest');
title('$C_{total}(t)/P_H(t)$');
grid on

Ph = SH+EH+AH+DH;
Nh = trapz(Ph,1)*da;
subplot(2,2,2), plot(t/365,(trapz(Ctot,1)*da)./Nh);
% axis_years(gca,tfinal)
title('$\int C_{total}(\alpha,t)d\alpha / N_H(t)$');
xlabel('time');
grid on

subplot(2,2,3), imagesc(t/365,a/365,Ctot./PH);
set(gca,'YDir','normal');
colorbar;
ylabel('age');
xlabel('time');
title('$C_{total}(\alpha,t)/P_H(\alpha,t)$');

subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/4)));
hold on;
subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/2)));
subplot(2,2,4), plot(a/365,Ctot(:,floor(3*nt/4)));
subplot(2,2,4), plot(a/365,Ctot(:,end));
xlabel('age')
legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
    ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)]);
title('$C_{total}(t)$');

%% Immunity breakdown
figure_setups; 
plot(a/365,Cac(:,end));
hold on;
plot(a/365,Cm(:,end));
xlabel('age (years)')
legend('Acquired','Maternal');
title('Immun dist.');
grid on
%% Impact of immunity on the sigmoids (rho, psi, phi)
% figure_setups;
% subplot(2,2,1), imagesc(t/365,a/365,sigmoid_prob(Ctot, 'phi'));
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('time');
% title('$\phi$');
% 
% subplot(2,2,2), imagesc(t/365,a/365,sigmoid_prob(Ctot, 'rho'));
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('time');
% title('$\rho$');
% 
% subplot(2,2,3), imagesc(t/365,a/365,sigmoid_prob(Ctot, 'psi'));
% set(gca,'YDir','normal');
% colorbar;
% ylabel('age');
% xlabel('time');
% title('$\psi$');
% 
% figure_setups;
% plot(a,P.rho)
% axis_years(gca,tfinal)
% title('$\rho$(age) at tfinal')
% grid on

%% Write output for convergence check
% Ctot_t = trapz(Ctot,1)*da; % total immunity in time
% Ctot_end = Ctot(:,end); % final distribution of immunity
% save(['Results/solution_',num2str(dt),'.mat'],'P','t','Ctot_t','Ctot_end')

%% Mosquito infection dynamics
figure_setups;
plot(t,SM,'b-'); hold on;
plot(t,EM,'-','Color',colour_r1);
plot(t,IM,'r-.');
plot(t,SM+EM+IM,'-.')
legend('SM','EM','IM','$N_M$');
title('mosquito population size by stages')
%axis_years(gca,tfinal); % change to x-axis to years if needed
grid on
% axis([0 tfinal 0 5])

%% final age distributions
% figure_setups;
% plot(a/365,SH(:,end),'-','Color',colour_mat1); hold on;
% plot(a/365,EH(:,end),'--','Color',colour_mat3);
% plot(a/365,AH(:,end),'-.','Color',colour_mat2);
% plot(a/365,DH(:,end),'-','Color',colour_mat7);
% legend('SH-age','EH-age','AH-age', 'DH-age');
% title('Age distributions by class (final time)')
% grid on
% 
%% Compare initial and final age distribution
% figure_setups;
% Ph = SH+EH+AH+DH;
% plot(a/365, Ph(:,end),'g-'); hold on;
% %plot(a/365,Ph(:,ceil(end/2)),'-','Color',colour_r1); 
% plot(a/365, P.n_tilde,'-','Color',colour_r2);
% plot(a/365, P.n_tilde*exp(P.p_hat*tfinal) ,'-.');
% legend('Final Age dist. (sim.)','Stable Age Dist. (initial)','Final Age dist. (theory)');
% title('Population Age Distributions')

%% R0 calculation
C_star = P.Lambda*(P.bm*P.bm*P.bh*P.bh*NM*NM)*P.betaM*P.sigma./...
    ((P.bm*NM*P.bh).^2).*(P.sigma+P.muM).*P.muM;
% NB the immunity functions are equal to 1/2 currently so the argument doesn't
% matter for them; functions only take scalar input right now
F_P = @(p,a) P.rho(1)*P.h*exp(-(p+P.rD)*a).*( P.h -P.h*exp((p+P.rD)*a)...
    -(p+P.rD).*exp((P.rD-P.h)*a) + p + P.rD*exp((p+P.rD)*a) )./((P.h+p)*(P.rD-P.h)*(p+P.rD));
% numerical overflow in F_P due to the term exp(age_max*P.h) appearing
G_P = @(p,a) (((1-P.rho(1))*P.h)./(P.h+p))*(1-exp(-(P.h+p)*a))+(1-P.phi(1))*P.rD*F_P(p,a);

H_P = @(p,a) exp(-(P.rA+p)*a).*da.*trapz(G_P(p,0:da:a).*exp((0:da:a).*(P.rA+p)));

zeta_P = @(p) C_star*trapz(exp(-P.muH_int).*(P.betaD.*F_P(p,age_max) + P.betaA.*H_P(p,age_max))) - 1;

% zeta_P(-0.1)
% zeta_P(0)
% zeta_P(0.1)













