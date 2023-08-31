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
P.vp0 = 0.8;
vp = zeros(size(P.a));
vp(vage_ind1:vage_ind2) = P.vp0;

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