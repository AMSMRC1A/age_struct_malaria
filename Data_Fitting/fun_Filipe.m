%% objective function for model calibration
function err = fun_Filipe(x)
global F
global P
% F(age(years), EIR) 
% =  immunity (~ Ctot) from RB's paper, Fig 5
% (~ rho) from Filipe's paper

P.phi_s_2 = x(1);
P.phi_r_2 = x(2);
P.rho_s_2 = x(3);
P.rho_r_2 = x(4); 
P.psi_s_2 = x(3);
P.psi_r_2 = x(4);

Malaria_parameters_transform;
nsamp = 10;
[~,ind1] = min(abs(P.a-0.5*365)); % start from 0.5 years old
[~,ind2] = min(abs(P.a-10*365)); % end at 10 years old
ind = round(linspace(ind1,ind2,nsamp)');

P.betaM = 0.2; % low EIR region ~ 25
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE','fsolve');
EIR1 = fit_EIR(SH,EH,DH,AH);
if EIR1<5; keyboard; end
x1 = P.a(ind)/365;
y1 = EIR1*ones(size(x1)); % aEIR
z1 = Ctot(ind)./P.PH_stable(ind); % final Ctot at EE

P.betaM = 0.5; % mid EIR ~ 100
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE','fsolve');
EIR2 = fit_EIR(SH,EH,DH,AH); 
x2 = P.a(ind)/365;
y2 = EIR2*ones(size(x2)); % aEIR
z2 = Ctot(ind)./P.PH_stable(ind); % final Ctot at EE

P.betaM = 1; % high EIR ~ 150
Malaria_parameters_transform;
[SH,EH,DH,AH,~,~,Ctot] = steady_state('EE','fsolve');
EIR3 = fit_EIR(SH,EH,DH,AH); 
x3 = P.a(ind)/365;
y3 = EIR3*ones(size(x3)); % aEIR
z3 = Ctot(ind)./P.PH_stable(ind); % final Ctot at EE

x = [x1;x2;x3];
y = [y1;y2;y3];
z = [z1;z2;z3]; 
z_data = F(x,y); % rho from data

z_samp = sigmoid_prob(z, 'rho'); % rho from samples

res = abs(z_data-z_samp);
w = ones(size(z_data));%./(res+eps);
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