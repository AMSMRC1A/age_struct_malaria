clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic

%% numerical config
P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero
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

P.np = 2; % number of patches

% model parameters - rates are in 1/day
baseline_Malaria_parameters_patch;

% allocation
% SH, EH, etc.. = cell averages
SH = NaN(na,nt,np); EH = NaN(na,nt,np); DH = NaN(na,nt,np); AH = NaN(na,nt,np);
SM = NaN(nt,np); EM = NaN(nt,np); IM = NaN(nt,np);

% per-person immunity level or \tilde(Ctotal) \tilde(Cm) \tilde(Cac) on overleaf
Ctot = NaN(na,nt,np); Cm = NaN(na,nt,np); Cac = NaN(na,nt,np);
NM = P.gM/P.muM*ones(nt,np);
NH = ones(nt,np);

% initial condition
% NH(1,np) NM(1,np) SH(na,np) SM(1,np) Cm(na,nt,np)
[SH(:,1,:),EH(:,1,:),DH(:,1,:),AH(:,1,:),SM(1,:),EM(1,:),IM(1,:)] = Malaria_IC(NH(1,:),NM(1,:));
[Cm(:,1,:),Cac(:,1,:),Ctot(:,1,:)] = Immunity_IC; % define initial immunity and related probability

%% time evolution
for n = 1:nt-1
    if mod(n,(nt-1)/10)==0
        disp(['progress = ',num2str(n/(nt-1)*100),'%']);
    end
    
    % human equations
    PH = SH(:,n,:)+EH(:,n,:)+DH(:,n,:)+AH(:,n,:); % PH(na,np) total human at age a, t = n, PH(na,np)
    NH(n,:) = trapz(PH)*da; % NH(nt,np) total human population at t=n; % trapz integrates over each column and returns a row vector
    NM(n,:) = SM(n,:)+EM(n,:)+IM(n,:); % SM(nt,np)  NM(nt,np)
    NH_star = NH(n,:)*P.pmat'; % NH_star(1,np)
    lamH = FOI_H(IM(n,:),NM(n,:),NH_star);  % lamH(1,np) force of infection for patch i
    lamH_star = lamH*P.pmat; % lamH_star(1,np), total FOI for people live in patch i
    % birth terms
    SH(1,n+1,:) = trapz(P.gH.*PH)*da; % PH(na,np) -> (1,np)
    EH(1,n+1,:) = 0;
    DH(1,n+1,:) = 0;
    AH(1,n+1,:) = 0;
    
    for ip = 1:P.np
        SH(2:end,n+1,ip) = (SH(1:end-1,n,ip)+P.dt*(P.phi(1:end-1,ip)*P.rD.*DH(1:end-1,n,ip)+P.rA*AH(1:end-1,n,ip))) ./ (1+(lamH_star(1,ip)+P.muH(2:end))*P.dt);
        EH(2:end,n+1,ip) = (EH(1:end-1,n,ip)+P.dt*lamH_star(1,ip)*SH(2:end,n+1,ip)) ./ (1+(P.h+P.muH(2:end))*P.dt);
        AH(2:end,n+1,ip) = ((1-P.dt*P.rA)*AH(1:end-1,n,ip)+P.dt*((1-P.rho(1:end-1,ip))*P.h.*EH(2:end,n+1,ip)+(1-P.phi(1:end-1,ip)).*P.rD.*DH(1:end-1,n,ip)))...
            ./ (1+P.dt*(P.psi(1:end-1,ip).*lamH_star(1,ip)+P.muH(2:end)));
        DH(2:end,n+1,ip) = ((1-P.dt*P.rD)*DH(1:end-1,n,ip)+P.dt*(P.rho(1:end-1,ip)*P.h.*EH(2:end,n+1,ip)+P.psi(1:end-1,ip).*lamH_star(1,ip).*AH(2:end,n+1,ip)))...
            ./ (1+P.dt*(P.muH(2:end)+P.muD(2:end)));
    end
    
    % mosquito equations
    [SM(n+1,:),EM(n+1,:),IM(n+1,:)] = mosquito_ODE(DH(:,n,:),AH(:,n,:),NH(n,:),NM(n,:));
    
    %% immunity equations (need to be modified for patch version)
    %     Cm(1,n+1) = P.m*trapz(P.gH.*PH.*Cac(:,n))*da/SH(1,n+1);
    %     Cac(1,n+1) = 0;
    %     % maternal immunity
    %     %n0 = min(n,na-1); % comment this formula for now, use implicit scheme
    %     %Cm(2:n0+1,n+1) = (Cm(1,1:n0))'.*exp(-a(2:n0+1)/P.dm); % k=1:n0
    %     %Cm(n0+2:end,n+1) = Cm(2:end-n0,1).*exp(-t(n+1)/P.dm);  % k=n0+1:end-1
    %     % acquired immunity - use Qn+1
    %     PHp1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1); % total human at age a, t = n+1
    %     NHp1 = trapz(PHp1)*da; % total human population at t=n;
    %     NMp1 = SM(1,n+1)+EM(1,n+1)+IM(1,n+1);
    %     [bHp1,~] = biting_rate(NHp1,NMp1);
    %     lamHp1 = FOI_H(bHp1,IM(1,n+1),NMp1);
    %     % neglecting disease induced mortality in Cac
    %     Qnp1 = f(lamHp1).*(P.cS*SH(2:end,n+1) + P.cE*EH(2:end,n+1) + P.cA*AH(2:end,n+1) ...
    %         + P.cD*DH(2:end,n+1)) + P.cV*P.v(2:end).*SH(2:end,n+1) ;
    %     Cac(2:end,n+1) = (Cac(1:end-1,n)+P.dt*Qnp1)./(1 + P.dt*(1./P.dac + P.muH(2:end)));
    %     Cm(2:end,n+1) = Cm(1:end-1,n)/(1+P.dt/P.dac);
    %     % Cm is now per person but Cac is pooled so need to multiply Cm by PH
    %     % to get total immunity contribution
    %     Ctot(:,n+1) = P.c1*Cac(:,n+1)+P.c2*Cm(:,n+1).*PHp1; % total immunity from acquired and maternal sources
    
    % update progression probability based on immunity Ctot
    Ctot(:,n+1,:) = Ctot(:,n,:);
    P.phi = sigmoid_prob(Ctot(:,n+1,:), 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Ctot(:,n+1,:), 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Ctot(:,n+1,:), 'psi'); % prob. AH -> DH
end
PH_final = reshape(SH(:,end,:)+EH(:,end,:)+DH(:,end,:)+AH(:,end,:),na,np); % PH_final(na,np) total human at age a, t = n
NH(end,:) = trapz(PH_final)*da;
toc
%% Population size versus time
figure_setups;
for ip = 1:np
    subplot(1,np,ip)
    Nh = (trapz(SH(:,:,ip),1)+trapz(EH(:,:,ip),1)+trapz(AH(:,:,ip),1)+trapz(DH(:,:,ip),1))*da;
    plot(t,trapz(SH(:,:,ip),1)*da,'-','Color',colour_mat1); hold on;
    plot(t,trapz(EH(:,:,ip),1)*da,'--','Color',colour_mat3);
    plot(t,trapz(AH(:,:,ip),1)*da,'-.','Color',colour_mat2);
    plot(t,trapz(DH(:,:,ip),1)*da,'-','Color',colour_mat7);
    plot(t,Nh,'-.k')
    legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$','Location','NorthWest');
    title(['Population size vs time (patch ',num2str(ip),')']);
    axis_years(gca,tfinal); % change to x-axis to years if needed
    grid on
    axis([0 tfinal 0 max(Nh)+0.1]);
end

%% Age profiles at tfinal
figure_setups;
for ip = 1:np
    subplot(1,np,ip)
    plot(a,SH(:,end,ip),'-','Color',colour_mat1); hold on;
    plot(a,EH(:,end,ip),'--','Color',colour_mat3);
    plot(a,AH(:,end,ip),'-.','Color',colour_mat2);
    plot(a,DH(:,end,ip),'-','Color',colour_mat7);
    plot(a,(SH(:,end,ip)+AH(:,end,ip)+EH(:,end,ip)+DH(:,end,ip)),'-.k');
    legend('SH','EH','AH', 'DH','PH');
    title(['Final Age Dist. (patch ',num2str(ip),')']);
    axis_years(gca,age_max); % change to x-axis to years if needed
    grid on
    axis([0 age_max 0 max(max(SH(:,end,ip)+AH(:,end,ip)+EH(:,end,ip)+DH(:,end,ip)))]);
end
%% Age proportions at tfinal
figure_setups;
for ip = 1:np
    subplot(1,np,ip)
    plot(a,SH(:,end,ip)./PH_final,'-','Color',colour_mat1); hold on;
    plot(a,EH(:,end,ip)./PH_final,'--','Color',colour_mat3);
    plot(a,AH(:,end,ip)./PH_final,'-.','Color',colour_mat2);
    plot(a,DH(:,end,ip)./PH_final,'-','Color',colour_mat7);
    plot(a,(SH(:,end,ip)+AH(:,end,ip)+EH(:,end,ip)+DH(:,end,ip))./PH_final(:,ip),'-.k');
    legend('SH','EH','AH', 'DH','PH');
    title(['Final Age Dist. Proportions (patch ',num2str(ip),')']);
    xlabel('age');
    axis_years(gca,age_max); % change to x-axis to years if needed
    grid on
    axis([0 age_max 0 1.1]);
end

%% Population proportions versus time
figure_setups;
for ip = 1:np
    subplot(1,np,ip)
    Nh = (trapz(SH(:,:,ip),1)+trapz(EH(:,:,ip),1)+trapz(AH(:,:,ip),1)+trapz(DH(:,:,ip),1))*da;
    plot(t,trapz(SH(:,:,ip),1)*da./Nh,'-','Color',colour_mat1); hold on;
    plot(t,trapz(EH(:,:,ip),1)*da./Nh,'--','Color',colour_mat3);
    plot(t,trapz(AH(:,:,ip),1)*da./Nh,'-.','Color',colour_mat2);
    plot(t,trapz(DH(:,:,ip),1)*da./Nh,'-','Color',colour_mat7);
    plot(t,(trapz(SH(:,:,ip),1)+trapz(EH(:,:,ip),1)+trapz(AH(:,:,ip),1)+trapz(DH(:,:,ip),1))*da./Nh,'-.k');
    legend('SH-age','EH-age','AH-age', 'DH-age','$N_H$');
    title(['Population proportions vs time (patch ',num2str(ip),')']);
    axis_years(gca,tfinal); % change to x-axis to years if needed
    xlabel('time');
    grid on
    axis([0 tfinal 0 1.1]);
end

%% Mosquito infection dynamics
figure_setups;
for ip = 1:np
    subplot(1,np,ip)
    plot(t,SM(:,ip),'b-'); hold on;
    plot(t,EM(:,ip),'-','Color',colour_r1);
    plot(t,IM(:,ip),'r-.');
    plot(t,SM(:,ip)+EM(:,ip)+IM(:,ip),'-.')
    legend('SM','EM','IM','$N_M$');
    title('mosquito population size by stages in patch i')
    %axis_years(gca,tfinal); % change to x-axis to years if needed
    grid on
    % axis([0 tfinal 0 5])
end
%% Immunity related figures
% for ip = 1:np
%     figure_setups;
%     PH = SH(:,:,ip)+EH(:,:,ip)+DH(:,:,ip)+AH(:,:,ip);
%     subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/4),ip)./PH(:,floor(nt/4)));
%     hold on;
%     subplot(2,2,1), plot(a/365,Ctot(:,floor(nt/2),ip)./PH(:,floor(nt/2)));
%     subplot(2,2,1), plot(a/365,Ctot(:,floor(3*nt/4),ip)./PH(:,floor(3*nt/4)));
%     subplot(2,2,1), plot(a/365,Ctot(:,end,ip)./PH(:,end),'-.');
%     % C_max = P.cS*SH(:,end,ip)/PH(:,end) + P.cE*SH(:,end,ip)/PH(:,end) + P.cA + P.cD;
%     % Im_bound = (P.dac/10).*(1 - exp(-a/P.dac));
%     % subplot(2,2,1), plot(a/365,Im_bound,'-.g');
%     xlabel('age')
%     legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
%         ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)],'Location','NorthWest');
%     title('$C_{total}(t)/P_H(t)$ in (patch i)');
%     grid on
%     
%     Nh = trapz(PH,1)*da;
%     subplot(2,2,2), plot(t/365,(trapz(Ctot(:,:,ip),1)*da)./Nh);
%     % axis_years(gca,tfinal)
%     title('$\int C_{total}(\alpha,t)d\alpha / N_H(t)$');
%     xlabel('time');
%     grid on
%     
%     subplot(2,2,3), imagesc(t/365,a/365,Ctot(:,:,ip)./PH);
%     set(gca,'YDir','normal');
%     colorbar;
%     ylabel('age');
%     xlabel('time');
%     title('$C_{total}(\alpha,t)/P_H(\alpha,t)$');
%     
%     subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/4),ip));
%     hold on;
%     subplot(2,2,4), plot(a/365,Ctot(:,floor(nt/2),ip));
%     subplot(2,2,4), plot(a/365,Ctot(:,floor(3*nt/4),ip));
%     subplot(2,2,4), plot(a/365,Ctot(:,end,ip),'-.');
%     xlabel('age')
%     legend(['t = ',num2str(tfinal/(4*365))],['t = ',num2str(tfinal/(2*365))],...
%         ['t = ',num2str(3*tfinal/(4*365))],['t = ',num2str(tfinal/365)]);
%     title('$C_{total}(t)$');
% end
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
% figure_setups;
% for ip = 1:np
%     subplot(1,np,ip)
%     plot(a,Cac(:,end,ip),'-.r');
%     hold on;
%     plot(a,Cm(:,end,ip),'-.b');
%     plot(a,P.c2*Cm(:,end,ip)+P.c1*Cac(:,end,ip),'-.k');
%     xlabel('age (years)')
%     legend('Acquired','Maternal','Total');
%     title('Immun dist. in patch i');
%     axis([0 age_max 0 max(max(Cm(:,end,ip)),max(Cac(:,end,ip)))*1.1]);
%     grid on
% end
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

%% Stability of DFE, only valid when q = 0
% if P.balance_fertility == 1
%     C_star = P.Lambda*(P.bm*P.bm*P.bh*P.bh*NM)*P.betaM*P.sigma./...
%         (((P.bm*NM + P.bh).^2).*(P.sigma+P.muM).*P.muM);
%     % NB the immunity functions are fixed equal to 1/2 currently so the argument doesn't
%     % matter for them; functions only take scalar input right now
%     F_P = @(p,a) P.rho(1)*P.h*( P.h*exp(-(p+P.rD).*a) - P.h...
%         -(p+P.rD).*exp(-(p+P.h).*a) + p*exp(-(p+P.rD).*a) + P.rD )./((P.h+p)*(P.rD-P.h)*(p+P.rD));
%     G_P = @(p,a) (((1-P.rho(1))*P.h)./(P.h+p)).*(1-exp(-(P.h+p).*a)) + (1-P.phi(1))*P.rD*F_P(p,a);
%     H_P = @(p,a) exp(-(P.rA+p)*a).*da.*trapz(G_P(p,0:da:a).*exp((0:da:a).*(P.rA+p)));
%     zeta_P = @(p) C_star*da*trapz(exp(-P.muH_int).*(P.betaD.*F_P(p,0:da:age_max)' + P.betaA.*H_P(p,0:da:age_max)')) - 1;
%     options = optimset('TolX',1e-12);
%     p0 = [-0.02 1]; % [LP,RP] -> need to take care when selecting the endpoints as zeta_P may  blow-up for negative p
%     p_star = fzero(zeta_P,p0,options);
%     disp(['p* = ',num2str(p_star)]);
%     if p_star < 0
%         disp('DFE is stable');
%     else
%         disp('DFE is unstable');
%     end
% end
