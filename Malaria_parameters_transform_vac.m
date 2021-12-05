function Malaria_parameters_transform_vac

global P

%% boosting vaccination functions
vb_fun = @(age) P.vb0.*ones(size(age));
vb = vb_fun(P.a);
P.vb = vb;
P.vb_fun = vb_fun;

%% RTS,S (protection) vaccination functions
age_range = 30; % # of days to finish vaccination [9 month, 9 month + age_range]
[~,vage_ind1] = min(abs(P.a-P.vage));
[~,vage_ind2] = min(abs(P.a-(P.vage+age_range)));
% NH = trapz(P.PH_stable)*P.da;
% [SH_EE,~,~,~,~,~,~] = steady_state('EE','numerical');
% if vage_ind1 == vage_ind2
%     P.vp0 = P.vyear/365/(P.NN*SH_EE(vage_ind1)/NH); % per-capita daily vaccine rate
% else
%     P.vp0 = P.vyear/365/(P.NN*trapz(SH_EE(vage_ind1:vage_ind2)*P.da)/NH); % per-capita daily vaccine rate
% end
P.vp0 = 0.8;
vp = zeros(size(P.a));
vp(vage_ind1:vage_ind2) = P.vp0;
% vp = zeros(size(P.a));
% vp = P.vyear/365*exppdf(P.a-9*30,15);
% vp_fun = @(age) P.vp0.*(age>=P.a(vage_ind1)).*(age<=P.a(vage_ind2));
%% approximation of theta at DFE - needed for analytical purpose: DFE, R0, bifurcation
% pi_fun = @(x) P.w+P.e.*vp_fun(x);
% pi_int_a = intf(pi_fun,P.a);
% pi_int_fun = @(x) interp1(P.a,pi_int_a,x);
% exp_pi_int_a = intf(@(x) P.w.*exp(pi_int_fun(x)),P.a);
% exp_pi_int = @(x) interp1(P.a,exp_pi_int_a,x);
% theta_fun = @(x) exp(-pi_int_fun(x)).*(1+exp_pi_int(x));
% theta0 = theta_fun(P.a);

%%

P.vp = vp;
% P.vp_fun = vp_fun;
% P.theta_fun = theta_fun;
% P.theta = theta0;


end

function intfx = intf(f,xs)
ns = length(xs);
fx = f(xs);
e = ones(ns,ns+1)/2;
e(2:end,1:end-1)=1;e(1,end)=0;
A = spdiags(e,-ns:0,ns,ns);
intfx= A*fx;
end