clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic

%% Numerical configuration and output settings
P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero
P.calc_R0 = 1; % 0 to ignore R0 calculation, 1 to compute it based on constant sigmoids
P.numerical_stability = 0; % 0 to numerical DFE stability calculation

tfinal = 200*365; % final time in days
age_max = 60*365; % max ages in days
P.age_max = age_max;
dt = 10; % time/age step size in days
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
P.halfmaximal = 1;

% model parameters - rates are in 1/day
baseline_Malaria_parameters;

% allocation
% SH, EH, etc.. = cell averages
SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt);
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt);

% per-person immunity level or \tilde(Ctotal) \tilde(Cm) \tilde(Cac) on Overleaf
Ctot = NaN(na,nt); Cm = NaN(na,nt); Cac = NaN(na,nt);
NM = P.gM/P.muM;
NH = ones(1,length(t));

% initial condition
[SH(:,1),EH(:,1),DH(:,1),AH(:,1),SM(1,1),EM(1,1),IM(1,1)] = Malaria_IC(NH(1),NM);
[Cm(:,1),Cac(:,1),Ctot(:,1)] = Immunity_IC; % initial immunity and related probability

%% time evolution
for n = 1:nt-1
    if mod(n,(nt-1)/5)==0
        disp(['progress = ',num2str(n/(nt-1)*100),'%']);
    end
    PH = SH(:,n)+EH(:,n)+DH(:,n)+AH(:,n); % total human at age a, t = n
    NH(n) = trapz(PH)*da; % total human population at t=n;
    NM = SM(1,n)+EM(1,n)+IM(1,n);
    [bH,~] = biting_rate(NH(n),NM);
    lamH = FOI_H(bH,IM(1,n),NM);  % force of infection at t=n
    
    % update progression probability based on immunity PER PERSON
    P.phi = sigmoid_prob(Ctot(:,n)./PH, 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Ctot(:,n)./PH, 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Ctot(:,n)./PH, 'psi'); % prob. AH -> DH
    
    % human birth terms
    SH(1,n+1) = trapz(P.gH.*PH)*da;
    EH(1,n+1) = 0;
    DH(1,n+1) = 0;
    AH(1,n+1) = 0;
    
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
    %n0 = min(n,na-1); % comment this formula for now, use implicit scheme
    %Cm(2:n0+1,n+1) = (Cm(1,1:n0))'.*exp(-a(2:n0+1)/P.dm); % k=1:n0
    %Cm(n0+2:end,n+1) = Cm(2:end-n0,1).*exp(-t(n+1)/P.dm);  % k=n0+1:end-1
    % acquired immunity - use Qn+1
    PHp1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1); % total human at age a, t = n+1
    NHp1 = trapz(PHp1)*da; % total human population at t=n;
    NMp1 = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    [bHp1,~] = biting_rate(NHp1,NMp1);
    lamHp1 = FOI_H(bHp1,IM(1,n+1),NMp1);
    % neglecting disease induced mortality in Cac
    Qnp1 = f(lamHp1).*(P.cS*SH(2:end,n+1) + P.cE*EH(2:end,n+1) + P.cA*AH(2:end,n+1) ...
        + P.cD*DH(2:end,n+1)) + P.cV*P.v(2:end).*SH(2:end,n+1) ;
    Cac(2:end,n+1) = (Cac(1:end-1,n)+P.dt*Qnp1)./(1 + P.dt*(1./P.dac + P.muH(2:end)));
    Cm(2:end,n+1) = Cm(1:end-1,n)/(1+P.dt/P.dac);
    % Cm is now per person but Cac is pooled so need to multiply Cm by PH
    % to get total immunity contribution
    Ctot(:,n+1) = P.c1*Cac(:,n+1)+P.c2*Cm(:,n+1).*PHp1; % total immunity from acquired and maternal sources
end
PH_final = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end); % total human at age a, t = n
NH(end) = trapz(PH_final)*da;
%% Population size versus time
figure_setups;
Nh = (trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da;
plot(t,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
plot(t,trapz(EH,1)*da,'--','Color',colour_mat3);
plot(t,trapz(AH,1)*da,'-.','Color',colour_mat2);
plot(t,trapz(DH,1)*da,'-','Color',colour_mat7);
plot(t,Nh,'-.k')
legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$','Location','NorthWest');
title('Population size vs time');
axis_years(gca,tfinal); % change to x-axis to years if needed
grid on
axis([0 tfinal 0 max(Nh)+0.1]);
%% Age profiles at tfinal
figure_setups;
plot(a,SH(:,end),'-','Color',colour_mat1); hold on;
plot(a,EH(:,end),'--','Color',colour_mat3);
plot(a,AH(:,end),'-.','Color',colour_mat2);
plot(a,DH(:,end),'-','Color',colour_mat7);
plot(a,(SH(:,end)+AH(:,end)+EH(:,end)+DH(:,end)),'-.k');
legend('SH','EH','AH', 'DH','PH');
title('Final Age Dist.');
axis_years(gca,age_max); % change to x-axis to years if needed
grid on
axis([0 age_max 0 max(max(SH(:,end)+AH(:,end)+EH(:,end)+DH(:,end)))]);
%% Age proportions at tfinal
figure_setups;
plot(a,SH(:,end)./PH_final,'-','Color',colour_mat1); hold on;
plot(a,EH(:,end)./PH_final,'--','Color',colour_mat3);
plot(a,AH(:,end)./PH_final,'-.','Color',colour_mat2);
plot(a,DH(:,end)./PH_final,'-','Color',colour_mat7);
plot(a,(SH(:,end)+AH(:,end)+EH(:,end)+DH(:,end))./PH_final,'-.k');
legend('SH','EH','AH', 'DH','PH');
title('Final Age Dist. Proportions');
xlabel('age');
axis_years(gca,age_max); % change to x-axis to years if needed
grid on
axis([0 age_max 0 1.1]);
%% Population proportions versus time
figure_setups;
plot(t,trapz(SH,1)*da./Nh,'-','Color',colour_mat1); hold on;
plot(t,trapz(EH,1)*da./Nh,'--','Color',colour_mat3);
plot(t,trapz(AH,1)*da./Nh,'-.','Color',colour_mat2);
plot(t,trapz(DH,1)*da./Nh,'-','Color',colour_mat7);
plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da./Nh,'-.k');
legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$');
title('Population proportions vs time');
axis_years(gca,tfinal); % change to x-axis to years if needed
xlabel('time');
grid on
axis([0 tfinal 0 1.1]);
%% Immunity related figures
figure_setups;
PH = SH+EH+DH+AH;
subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/4))./PH(:,floor(nt/4)));
hold on;
subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/2))./PH(:,floor(nt/2)));
subplot(2,2,1), plot(a/365,Ctot(:,floor(3*nt/4))./PH(:,floor(3*nt/4)));
subplot(2,2,1), plot(a/365,Ctot(:,end)./PH(:,end),'-.');
% C_max = P.cS*SH(:,end)/PH(:,end) + P.cE*SH(:,end)/PH(:,end) + P.cA + P.cD;
% Im_bound = (P.dac/10).*(1 - exp(-a/P.dac));
% subplot(2,2,1), plot(a/365,Im_bound,'-.g');
xlabel('age')
legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
    ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)],'Location','NorthWest');
title('$C_{total}(t)/P_H(t)$');
grid on

Nh = trapz(PH,1)*da;
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
subplot(2,2,4), plot(a/365,Ctot(:,end),'-.');
xlabel('age')
legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
    ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)]);
title('$C_{total}(t)$');
%% Immunity calculations at steady state
% C_SS = zeros(1,na);
% for i = 1:na
%     C_SS(i) = f(lamHp1)*da*exp(-a(i)/P.dac)*trapz(exp(a(1:i)./P.dac).*(P.cS*SH(1:i,end) + P.cE*EH(1:i,end)...
%         + P.cA*AH(1:i,end) + P.cD*DH(1:i,end))./PH_final(1:i));
% end
% figure_setups;
% plot(a,C_SS,'-.r');
% title('Steady state acquired immunity');
% axis_years(gca,age_max);
% xlabel('age')
% axis([0 age_max 0 max(C_SS)*1.1]);
% grid on;
%% Immunity breakdown
figure_setups;
plot(a,Cac(:,end),'-.r');
hold on;
plot(a,Cm(:,end),'-.b');
plot(a,P.c2*Cm(:,end)+P.c1*Cac(:,end),'-.k');
xlabel('age (years)')
legend('Acquired','Maternal','Total');
title('Immun dist.');
axis([0 age_max 0 max(max(Cm(:,end)),max(Cac(:,end)))*1.1]);
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
% figure_setups;
% plot(t,SM,'b-'); hold on;
% plot(t,EM,'-','Color',colour_r1);
% plot(t,IM,'r-.');
% plot(t,SM+EM+IM,'-.')
% legend('SM','EM','IM','$N_M$');
% title('mosquito population size by stages')
%axis_years(gca,tfinal); % change to x-axis to years if needed
%grid on
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

%% DFE stability, R0 calculation based on constant sigmoids
% Theoretical R0 calculation
if P.balance_fertility == 1 && P.calc_R0 == 1
    C_star = P.Lambda*(P.bm*P.bm*P.bh*P.bh*NM)*P.betaM*P.sigma./...
        (((P.bm*NM + P.bh).^2).*(P.sigma+P.muM).*P.muM);
    % NB the immunity functions are assumed constant in this calculation
    F_P = @(p,a) -P.rho(1)*P.h*( (p+P.h)*exp(-(p+P.rD).*a) - P.h...
        -(p+P.rD).*exp(-(p+P.h).*a) + P.rD )./((P.h+p)*(P.h-P.rD)*(p+P.rD));
    G_P = @(p,a) (((1-P.rho(1))*P.h)./(P.h+p)).*(1-exp(-(P.h+p).*a)) + (1-P.phi(1))*P.rD*F_P(p,a);
    H_P = @(p,a) (((1-P.rho(1))*P.h)./(P.h+p)).*( ( -P.h + (P.h+p)*exp(-(p+P.rA)*a) + P.rA - (p + P.rA)*exp(-(P.h+p)*a) )./((p + P.rA)*(P.rA - P.h)) )...
        + ((1-P.phi(1))*P.rD*P.rho(1)*P.h/((P.h+p)*(P.rD-P.h)*(p+P.rD)))*...
        ( ((P.h -P.rD)/(p+P.rA)-(p+P.rD)/(P.h-P.rA)+(P.h+p)/(P.rD-P.rA))*exp(-(p+P.rA)*a) + ...
        (p+P.h)*exp(-(p+P.rD)*a)/(P.rA-P.rD) + (P.rD-P.h)/(P.rA+p) + (p+P.rD)*exp(-(P.h+p)*a)/(P.h-P.rA) );
    zeta_P = @(p) C_star*da*trapz(exp(-P.muH_int).*(P.betaD.*F_P(p,0:da:age_max)' + P.betaA.*H_P(p,0:da:age_max)'));
    R0 = zeta_P(0);
    disp(['R0 = ',num2str(R0)]);
    if R0 < 1
        disp('R0 is less than 1; DFE is stable');
    else
        disp('R0 is greater than 1; DFE is unstable');
    end
end
%% Numerical DFE stability calculation, based on constant sigmoids
if P.balance_fertility == 1 && P.numerical_stability == 1
    % Compute the Jacobian and stability of the DFE numerically
    [bH,bM] = biting_rate(1,P.gM/P.muM);
    Lambda_M = @(x) bM*P.Lambda*da*trapz(exp(-P.muH_int).*(P.betaD*x(:,3)+P.betaA*x(:,4)));
    Lambda_H = @(x) bH*P.betaM*(P.sigma/(P.sigma+P.muM))*((Lambda_M(x))/(Lambda_M(x) + P.muM));
    F_PSI = @(x) [[-x(1,1)+1; -Lambda_H(x)*x(2:end,1) + P.phi(1)*P.rD*x(2:end,3) + P.rA*x(2:end,4) - diff(x(:,1))/da];...
        [-x(1,2); Lambda_H(x)*x(2:end,1) - P.h*x(2:end,2) - diff(x(:,2))/da];...
        [-x(1,3); P.rho(1)*P.h*x(2:end,2) + P.psi(1)*Lambda_H(x)*x(2:end,4) - P.rD*x(2:end,3) - diff(x(:,3))/da];...
        [-x(1,4); (1-P.rho(1))*P.h*x(2:end,2) - P.psi(1)*Lambda_H(x)*x(2:end,4) + (1-P.phi(1))*P.rD*x(2:end,3) - P.rA*x(2:end,4) - diff(x(:,4))/da]];
    error_tolerance = 1e-20;
    options = optimoptions('fsolve','Display','none','OptimalityTolerance',...
        error_tolerance);
    P1 = [ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1)];
    [eq_age,fval,exitflag,output,jacobian] = fsolve(F_PSI,P1,options);
    % get the eigenvalues of the Jacobian and plot the linearized spectrum
    temp_eig = sort(real(eig(jacobian)),'descend');
    re_max = max(temp_eig);
    if re_max < 0
        disp(['max real part of eigenvalues at DFE = ',num2str(re_max,'%10.6f'), '; DFE is stable']);
    else
        disp(['max real part of eigenvalues = ',num2str(re_max,'%10.6f'), '; DFE is unstable']);
    end
end
%%
toc