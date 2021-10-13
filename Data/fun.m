function err = fun(x)
global F
global P
% F(age(years), EIR) = susceptibility (~ rho) from Filipe's paper

P.c_rho = x(1);
P.c_psi = x(2);
P.c_phi = x(3);
P.k_rho = x(4);
P.k_psi = x(5);
P.k_phi = x(6);
Malaria_parameters_transform;

[~,ind1] = min(abs(P.a-0.1*365));
[~,ind2] = min(abs(P.a-1*365));
[~,ind3] = min(abs(P.a-5*365));
[~,ind4] = min(abs(P.a-10*365));
ind = round([linspace(ind1,ind2,5),linspace(ind2,ind3,5),linspace(ind3,ind4,5)]);
x = P.a(ind)/365;

P.betaM = 0.006; % low EIR region 
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
y1 = fit_EIR(SH,EH,DH,AH)*ones(size(x)); % aEIR
z1 = sigmoid_prob(Ctot(ind)./P.PH_stable(ind), 'rho'); % final rho function at EE

P.betaM = 0.02; % mid EIR 
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
y2 = fit_EIR(SH,EH,DH,AH)*ones(size(x)); % aEIR
z2 = sigmoid_prob(Ctot(ind)./P.PH_stable(ind), 'rho'); % final rho function at EE

P.betaM = 0.25; % high EIR 
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE');
y3 = fit_EIR(SH,EH,DH,AH)*ones(size(x)); % aEIR
z3 = sigmoid_prob(Ctot(ind)./P.PH_stable(ind), 'rho'); % final rho function at EE

x = [x;x;x];
y = [y1;y2;y3];
z = [z1;z2;z3];

res = abs(F(x,y)-z);
w = ones(size(z))./(res+eps);
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