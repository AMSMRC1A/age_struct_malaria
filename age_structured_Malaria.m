function [SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria()
global P

na = P.na;
nt = P.nt;
da = P.da;
t = P.t; 

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
    if mod(n,(nt-1)/5)==0
        disp(['progress = ',num2str(n/(nt-1)*100),'%']);
    end
    PH = SH(:,n)+EH(:,n)+DH(:,n)+AH(:,n); % total human at age a, t = n
    NH(n) = trapz(PH)*da; % total human population at t=n;
    NM = SM(1,n)+EM(1,n)+IM(1,n);
    [bH,~] = biting_rate(NH(n),NM);
    lamH = FOI_H(bH,IM(1,n),NM);  % force of infection at t=n
    
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
    % update progression probability based on immunity Ctot
    P.phi = sigmoid_prob(Ctot(:,n+1), 'phi'); % prob. of DH -> RH
    P.rho = sigmoid_prob(Ctot(:,n+1), 'rho'); % prob. of severely infected, EH -> DH
    P.psi = sigmoid_prob(Ctot(:,n+1), 'psi'); % prob. AH -> DH
end

end