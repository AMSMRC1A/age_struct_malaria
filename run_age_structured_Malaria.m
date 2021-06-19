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
P.t = t;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;

% model parameters - rates are in 1/day
baseline_Malaria_parameters;

NM = P.gM/P.muM;
NH = 1;

% initial condition 
[SH,EH,DH,AH,RH,SM,EM,IM] = Malaria_IC(NH,NM); 
[Cm,Cac,Cs] = Immunity_IC; % initial immunity and related probability
CmIC = Cm; % record Initial condition

% allocation
CmBC = NaN(1,nt); % record boundary condition
CmBC(1) = Cm(1);
rho_ave = NaN(1,nt);
SHp1 = NaN(size(SH)); EHp1 = SHp1; DHp1 = SHp1; AHp1 = SHp1; RHp1 = SHp1;
Cmp1 = NaN(size(Cm)); Cacp1 = Cmp1; Csp1 = Cmp1; 
rho_ave(1,1) = mean(P.rho);

% create file & save IC
FileSave(SH,EH,DH,AH,RH,SM,EM,IM,Cm,Cac,Cs,'w');

%% time evolution
for n = 1:nt-1
    if mod(n,(nt-1)/10)==0
        display(['progress = ',num2str(n/(nt-1)*100),'%'])
    end
    NHa = SH+EH+DH+AH+RH; % total human at age a, t = n
    NH = trapz(NHa)*da; % total human population at t=n;
    
    NM = SM+EM+IM;
    [bH,~] = biting_rate(NH,NM); 
    lamH = FOI_H(bH,IM,NM);  % force of infection at t=n
    
    % human birth terms
    SHp1(1) = trapz(P.gH.*NHa)*da;
    EHp1(1) = 0;
    DHp1(1) = 0;
    AHp1(1) = 0;
    RHp1(1) = 0;
    
    SHp1(2:end) = (SH(1:end-1)+P.dt*P.w(2:end).*RH(1:end-1)) ./ (1+(lamH+P.v(2:end)+P.muH(2:end))*P.dt);
    EHp1(2:end) = (EH(1:end-1)+P.dt*lamH*SHp1(2:end)) ./ (1+(P.h+P.muH(2:end))*P.dt);
    % solve DH and AH together
        num = (1+P.dt*(P.rA+ P.muH(2:end))).*(DH(1:end-1)+P.dt*P.h*P.rho(2:end).*EHp1(2:end))+...
            P.dt*(AH(1:end-1)+DH(1:end-1)+P.dt*P.h*EHp1(2:end)).*P.psi(2:end)*lamH;
        den = (1+P.dt*(P.rA+ P.muH(2:end))).*(1+P.dt*(P.rD+P.muD(2:end)+P.muH(2:end)))+...
            P.dt*(1+P.dt*(P.muD(2:end)+P.muH(2:end)+P.rD*P.phi(2:end))).*P.psi(2:end)*lamH;
    DHp1(2:end) = num./den;
	AHp1(2:end) = (AH(1:end-1)+P.dt*((1-P.rho(2:end))*P.h.*EHp1(2:end)+(1-P.phi(2:end))*P.rD.*DHp1(2:end))) ./ (1+P.dt*(P.psi(2:end)*lamH+P.rA+P.muH(2:end)));
	RHp1(2:end) = ((1-P.w(2:end)*P.dt).*RH(1:end-1)+P.dt*(P.phi(2:end)*P.rD.*DHp1(2:end)+P.rA*AHp1(2:end)+P.v(2:end).*SHp1(2:end))) ./ (1+P.dt*P.muH(2:end));
        
    % adjust mosquito infection accordingly - use tn level!
    [SMp1,EMp1,IMp1] = mosquito_ODE(DH,AH,NH,NM);
    
    % immunity at age = 0
    Cmp1(1) = P.m*trapz(P.gH.*Cs.*NHa)/NH*da; % why is Cs in the integral here?
    Cacp1(1) = 0;
    CmBC(n+1) = Cmp1(1);
    % maternal immunity
    n0 = min(n,na-1);
    Cmp1(2:n0+1) = (CmBC(1:n0))'.*exp(-a(2:n0+1)/P.dm); % k=1:n0, use Boundary
    Cmp1(n0+2:end) = CmIC(2:end-n0).*exp(-t(n+1)/P.dm);  % k=n0+1:end-1, use Initial
    % acquired immunity
    NHap1 = SHp1+EHp1+DHp1+AHp1+RHp1; % total human at age a, t = n+1
    NHp1 = trapz(NHap1)*da; % total human population at t=n;
    NMp1 = SMp1+EMp1+IMp1;
    [bHp1,~] = biting_rate(NHp1,NMp1);
    lamHp1 = FOI_H(bHp1,IMp1,NMp1);   
    Qnp1 = lamHp1*(P.w1*SHp1(2:end) + P.w2*EHp1(2:end) + P.w3*AHp1(2:end) ...
        + P.w4*DHp1(2:end) + P.w5*RHp1(2:end))/NHp1;
    Cacp1(2:end) = (Cac(1:end-1)+P.dt*Qnp1)./(1+1/P.ds*P.dt); % use Qn+1
    Csp1 = P.c1*Cacp1+P.c2*Cmp1; % total immunity from acquired and maternal sources
    
    % update progression probability based on immunity Cs
    P.phi = sigmoid_prob(Csp1, 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Csp1, 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Csp1, 'psi'); % prob. AH -> DH
    rho_ave(1,n+1) = mean(P.rho);
    
    % file output
    FileSave(SHp1,EHp1,DHp1,AHp1,RHp1,SMp1,EMp1,IMp1,Cmp1,Cacp1,Csp1,'a')
    
    % update variables
    SH = SHp1; EH = EHp1; DH = DHp1; AH = AHp1; RH = RHp1; 
    SM = SMp1; EM = EMp1; IM = IMp1;
    Cm = Cmp1; Cac = Cacp1; Cs = Csp1;
end
save('Output/config.mat','P')
toc


