% close all
clear all
clc
format long
global P
global a 
global colour_r1 colour_r2

tic

% numerical config
tfinal = 400; % final time
nt = tfinal+1; % number of subdivisions in time/age
age_max = 50*365;
dt = tfinal/(nt-1); % grid size for time
da = dt; % grid size for age
t = 0:dt:tfinal; % time mesh
a = (0:da:age_max)'; % age mesh
na = length(a);

% model parameters
baseline_Malaria_parameters;

% allocation
SH = NaN(na,nt); EH = NaN(na,nt); DH = NaN(na,nt); AH = NaN(na,nt); RH = NaN(na,nt);
SM = NaN(1,nt); EM = NaN(1,nt); IM = NaN(1,nt); 
Cs = NaN(na,nt);
NM = P.gM/P.muM;
NH = 1;%trapz(P.gH./P.muH)*da; % need to parameterize the birth-death rates, so that the ratio M to H is reasonable

% initial condition 
[SH(:,1),EH(:,1),DH(:,1),AH(:,1),RH(:,1),SM(:,1),EM(:,1),IM(:,1)] = Malaria_IC(NH,NM); 
Cs(:,1) = Immunity_IC; % initial immunity and related probability
    
%% time evolution
for n = 1:nt-1
    NHa = SH(:,n)+EH(:,n)+DH(:,n)+AH(:,n)+RH(:,n); % total human at age a, t = n
    NH = sum(NHa); % total human population at t=n;
    NM = sum(SM(:,n)+EM(:,n)+IM(:,n));
    [bH,~] = biting_rate(NH,NM); 
    lamH = FOI_H(bH,IM(:,n),NM);  % force of infection at t=n
    
    % human birth terms
    SH(1,n+1) = da*trapz(P.gH.*NHa);
    EH(1,n+1) = 0;
    DH(1,n+1) = 0;
    AH(1,n+1) = 0;
    RH(1,n+1) = 0;
    % maternal immunity at age = 0
    Cs(1,n+1) = P.Cm0; 
    
    for k = 1:na-1
        SH(k+1,n+1) = (SH(k,n)+dt*P.w(k+1)*RH(k,n)) / (1+(lamH+P.v(k+1)+P.muH(k+1))*dt);
        EH(k+1,n+1) = (EH(k,n)+dt*lamH*SH(k+1,n+1)) / (1+(P.h+P.muH(k+1))*dt);
        % solve DH and AH together
        num = (1+dt*(P.rA+ P.muH(k+1)))*(DH(k,n)+dt*P.h*P.rho*EH(k+1,n+1))+...
            dt*(AH(k,n)+DH(k,n)+dt*P.h*EH(k+1,n+1))*P.psi*lamH;
        den = (1+dt*(P.rA+ P.muH(k+1)))*(1+dt*(P.rD+P.muD(k+1)+P.muH(k+1)))+...
            dt*(1+dt*(P.muD(k+1)+P.muH(k+1)+P.rD*P.phi))*P.psi*lamH;
        DH(k+1,n+1) = num/den;
        AH(k+1,n+1) = (AH(k,n)+dt*((1-P.rho)*P.h*EH(k+1,n+1)+(1-P.phi)*P.rD*DH(k+1,n+1))) / (1+dt*(P.psi*lamH+P.rA+P.muH(k+1)));
        RH(k+1,n+1) = (1-P.w(k+1)*dt)*RH(k,n)+dt*(P.phi*P.rD*DH(k+1,n+1)+P.rA*AH(k+1,n+1)+P.v(k+1)*SH(k+1,n+1))/(1+dt*P.muH(k+1));
        % constraint on timestep for positivity: (1-P.w(k+1)*dt) > 0               
              
        Qn = P.c1*lamH*SH(k,n)/NH + P.c2*lamH*AH(k,n)/NH;
        Cs(k+1,n+1) = (Cs(k,n)+dt*(Qn+(1/P.ds-1)*P.Cm(k+1)))/(1+1/P.ds*dt);       
    end
    
    % adjust mosquito infection accordingly
    NHap1 = SH(:,n+1)+EH(:,n+1)+DH(:,n+1)+AH(:,n+1)+RH(:,n+1); % total human at age a, t = n+1
    NHp1 = sum(NHap1); % total human population at t=n+1;
    NMp1 = sum(SM(:,n+1)+EM(:,n+1)+IM(:,n+1));
    [SM(:,n+1),EM(:,n+1),IM(:,n+1)] = mosquito_ODE(DH(:,n+1),AH(:,n+1),NHp1,NMp1);
    
    % update progression probability based on immunity Cs
    P.phi = sigmoid_prob(Cs(:,n+1), 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Cs(:,n+1), 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Cs(:,n+1), 'psi'); % prob. AH -> DH

end

%% Plotting
figure_setups;
plot(t,sum(SH,1),'b-'); hold on;
plot(t,sum(EH,1),'-','Color',colour_r1);
plot(t,sum(AH,1),'-','Color',colour_r2);
plot(t,sum(DH,1),'r-');
plot(t,sum(RH,1),'g-');
plot(t,sum(SH,1)+sum(EH,1)+sum(AH,1)+sum(DH,1)+sum(RH,1))
legend('SH-age','EH-age','AH-age', 'DH-age','RH-age','$N_H$');
title('human')
grid on
toc;

figure_setups;
plot(t,SM,'b-'); hold on;
plot(t,EM,'-','Color',colour_r1);
plot(t,IM,'r-');
plot(t,SM+EM+IM)
legend('SM','EM','IM','$N_M$');
title('mosquitoes')
grid on
toc;