clear all
clc
format long
global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

tic
%% numerical config
P.balance_fertility = 1; % 0 to keep original fertility, 1 for balanced birth rate so that pop. growth is zero
P.solution_plots = 0; % toggle solution plots
P.sigmoid_surfs = 0; % toggle surfaces from the sigmoids
P.age_aEIR_surf = 1;

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
baseline_Malaria_parameters; % set parameters

halfmaximal_var = [1];
gM_var = [0.01 0.05 0.1 0.5 1]; % NB muM = 1/10, NM = gM/muM baseline
final_immunity = zeros(na,length(gM_var));
final_EIR = zeros(1,length(gM_var));
for ii = 1:length(halfmaximal_var)
    P.halfmaximal = halfmaximal_var(ii);
    for jj = 1:length(gM_var)
        P.gM = gM_var(jj);
        %% This section is one simulation
        % allocation
        SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt);
        SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt);
        % per-person immunity level or \tilde(Ctotal) \tilde(Cm) \tilde(Cac) on overleaf
        Ctot = NaN(na,nt); Cm = NaN(na,nt); Cac = NaN(na,nt);
        NM = P.gM/P.muM;
        NH = ones(1,length(t));
        % initial condition
        [SH(:,1),EH(:,1),DH(:,1),AH(:,1),SM(1,1),EM(1,1),IM(1,1)] = Malaria_IC(NH(1),NM);
        [Cm(:,1),Cac(:,1),Ctot(:,1)] = Immunity_IC; % initial immunity and related probability
        P.gM = gM_var(jj); % need to reset this again after ICs get set
        % time evolution
        for n = 1:nt-1
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
                + P.cD*DH(2:end,n+1)) + P.cV*P.v(2:end).*SH(2:end,n+1);
            Cac(2:end,n+1) = (Cac(1:end-1,n)+P.dt*Qnp1)./(1 + P.dt*(1./P.dac + P.muH(2:end)));
            Cm(2:end,n+1) = Cm(1:end-1,n)/(1+P.dt/P.dac);
            % Cm is now per person but Cac is pooled so need to multiply Cm by PH
            % to get total maternal immunity contribution to the pool
            Ctot(:,n+1) = P.c1*Cac(:,n+1)+P.c2*Cm(:,n+1).*PHp1; % total immunity from acquired and maternal sources
        end
        %% Record data for the age vs aEIR surface to mirror the R-B paper
        if P.age_aEIR_surf == 1
            final_immunity(:,jj) = Ctot(:,end)./PHp1;
            [bH,~] = biting_rate(NH(end),NM);
            lamH = FOI_H(bH,IM(1,end),NM);
            final_EIR(1,jj) = lamH/P.betaM; 
            disp(['EIR = ',num2str(final_EIR(1,jj),'%10.6f')]);
        end
        %% Plot the sigmoids
        if P.sigmoid_surfs == 1
            figure_setups;
            % commented plots are for cross sections at the final time
            %     dI = max(Ctot(:,end)./PHp1)/1000;
            %     dom = 0:dI:max(Ctot(:,end)./PHp1);
            %     subplot(2,2,1), plot(dom,sigmoid_prob(dom, 'phi')); title('$\phi$'); xlabel('immunity pp'); ylabel('prob.');
            %     subplot(2,2,2), plot(dom,sigmoid_prob(dom, 'rho')); title('$\rho$'); xlabel('immunity pp'); ylabel('prob.');
            %     subplot(2,2,3), plot(dom,sigmoid_prob(dom, 'psi')); title('$\psi$'); xlabel('immunity pp'); ylabel('prob.');
            %     grid on
            PH = SH+EH+DH+AH;
            subplot(2,2,1), imagesc(t/365,a/365,sigmoid_prob(Ctot./PH, 'phi'),[0 1]); colorbar; title('$\phi$'); xlabel('time')
            set(gca,'YDir','normal');
            grid on
            subplot(2,2,2), imagesc(t/365,a/365,sigmoid_prob(Ctot./PH, 'rho'),[0 1]); colorbar; title('$\rho$'); xlabel('time')
            set(gca,'YDir','normal');
            grid on
            subplot(2,2,3), imagesc(t/365,a/365,sigmoid_prob(Ctot./PH, 'psi'),[0 1]); colorbar; title('$\psi$'); xlabel('time')
            set(gca,'YDir','normal');
            grid on
        end
        %% Plot the results
        if P.solution_plots == 1
            PH_final = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end); % total human at age a, t = n
            NH(end) = trapz(PH_final)*da;
            % Population size versus time
            figure_setups;
            Nh = (trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da;
            subplot(2,2,1), plot(t,trapz(SH,1)*da,'-','Color',colour_mat1); hold on;
            subplot(2,2,1), plot(t,trapz(EH,1)*da,'--','Color',colour_mat3);
            subplot(2,2,1), plot(t,trapz(AH,1)*da,'-.','Color',colour_mat2);
            subplot(2,2,1), plot(t,trapz(DH,1)*da,'-','Color',colour_mat7);
            subplot(2,2,1), plot(t,Nh,'-.k')
            legend('$S_H$','$E_H$','$A_H$', '$D_H$','$N_H$');
            title('Population size');
            axis_years(gca,tfinal); % change to x-axis to years if needed
            grid on
            axis([0 tfinal 0 max(Nh)+0.1]);
            % Age profiles at tfinal
            subplot(2,2,2), plot(a,SH(:,end),'-','Color',colour_mat1); hold on;
            subplot(2,2,2), plot(a,EH(:,end),'--','Color',colour_mat3);
            subplot(2,2,2), plot(a,AH(:,end),'-.','Color',colour_mat2);
            subplot(2,2,2), plot(a,DH(:,end),'-','Color',colour_mat7);
            subplot(2,2,2), plot(a,(SH(:,end)+AH(:,end)+EH(:,end)+DH(:,end)),'-.k');
            title('Final Age Dist.');
            axis_years(gca,age_max); % change to x-axis to years if needed
            grid on
            axis([0 age_max 0 max(max(SH(:,end)+AH(:,end)+EH(:,end)+DH(:,end)))]);
            % Population proportions versus time
            subplot(2,2,3), plot(t,trapz(SH,1)*da./Nh,'-','Color',colour_mat1); hold on;
            subplot(2,2,3), plot(t,trapz(EH,1)*da./Nh,'--','Color',colour_mat3);
            subplot(2,2,3), plot(t,trapz(AH,1)*da./Nh,'-.','Color',colour_mat2);
            subplot(2,2,3), plot(t,trapz(DH,1)*da./Nh,'-','Color',colour_mat7);
            subplot(2,2,3), plot(t,(trapz(SH,1)+trapz(EH,1)+trapz(AH,1)+trapz(DH,1))*da./Nh,'-.k');
            title('Proportions');
            axis_years(gca,tfinal); % change to x-axis to years if needed
            xlabel('time');
            grid on
            axis([0 tfinal 0 1.1]);
            % Immunity related figures
            PH = SH+EH+DH+AH;
            subplot(2,2,4), plot(a/365,Ctot(:,end)./PH(:,end),'-.');
            % C_max = P.cS*SH(:,end)/PH(:,end) + P.cE*SH(:,end)/PH(:,end) + P.cA + P.cD;
            % Im_bound = (P.dac/10).*(1 - exp(-a/P.dac));
            % subplot(2,2,1), plot(a/365,Im_bound,'-.g');
            xlabel('age')
            title('$C_{total}(t)/P_H(t)$');
            grid on
        end
    end
end
%% Plot the age vs aEIR surface to mirror the R-B paper
if P.age_aEIR_surf == 1
    imagesc(a/365,final_EIR,final_immunity'); title('Immunity levels');
    set(gca,'YDir','normal');
    grid on
    colorbar;
end
%%
toc