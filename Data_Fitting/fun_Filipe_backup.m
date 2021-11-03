%% objective function for model calibration %% backup, not used
function err = fun_Filipe(x)
global F
global P
% F(age(years), EIR) = susceptibility (~ rho) from Filipe's paper

P.phi_t_2 = x(1);
P.phi_s_2 = x(2);
P.rho_t_2 = x(3);
P.rho_s_2 = x(4); 
P.psi_t_2 = x(5);
P.psi_s_2 = x(6);

Malaria_parameters_transform;
nsamp = 10;
[~,ind1] = min(abs(P.a-0.1*365));
[~,ind2] = min(abs(P.a-1*365));
ind = round(linspace(ind1,ind2,5)');

P.betaM = 0.008; % low EIR region ~ 25
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
EIR1 = fit_EIR(SH,EH,DH,AH);
ind1 = [randsample(length(P.a),nsamp,true,F(P.a/365,EIR1*ones(size(P.a))));ind];
x1 = P.a(ind1)/365;
y1 = EIR1*ones(size(x1)); % aEIR
z1 = sigmoid_prob(Ctot(ind1)./P.PH_stable(ind1), 'rho'); % final rho function at EE

P.betaM = 0.02; % mid EIR ~ 100
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
EIR2 = fit_EIR(SH,EH,DH,AH);
ind2 = [randsample(length(P.a),nsamp,true,F(P.a/365,EIR2*ones(size(P.a))));ind];
x2 = P.a(ind2)/365;
y2 = EIR2*ones(size(x2)); % aEIR
z2 = sigmoid_prob(Ctot(ind2)./P.PH_stable(ind2), 'rho'); % final rho function at EE

P.betaM = 0.25; % high EIR ~ 150
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
EIR3 = fit_EIR(SH,EH,DH,AH);
ind3 = [randsample(length(P.a),nsamp,true,F(P.a/365,EIR3*ones(size(P.a))));ind];
x3 = P.a(ind3)/365;
y3 = EIR3*ones(size(x3)); % aEIR
z3 = sigmoid_prob(Ctot(ind3)./P.PH_stable(ind3), 'rho'); % final rho function at EE

x = [x1;x2;x3];
y = [y1;y2;y3];
z = [z1;z2;z3];

res = abs(F(x,y)-z);
w = ones(size(z));%./(res+eps);
err = sum(w.*(res.^2));

end

function EIR = fit_EIR(SH,EH,DH,AH)
global P
NH = trapz(SH+EH+DH+AH)*P.da;
NM = P.gM/P.muM;
[bH,bM] = biting_rate(NH,NM);
Lambda_M = bM*trapz(P.betaD*DH + P.betaA*AH)*P.da;
IM_frac_EE = P.sigma/(P.sigma+P.muM)*(Lambda_M/(Lambda_M + P.muM));
EIR = bH*IM_frac_EE*365; % annual EIR
end