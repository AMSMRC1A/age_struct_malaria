%% sigmoidal function for converting immunity to probabilties
% return function handle
function fun = sigmoid_prob_fun(lprob)
global P

sigmoid_inc = @(f_0, f_1, x, s_2, r_2) f_0 + (f_1-f_0)./(1 + exp(-(x-s_2)/r_2));
sigmoid_dec = @(f_0, f_1, x, s_2, r_2) f_1 + (f_0-f_1)./(1 + exp(-(x-s_2)/r_2));

% NB set f1 = f0 to get a constant function (equal to f0 everywhere)

switch lprob
    case 'phi'
        fun = @(x) sigmoid_inc(P.phi_f_0, P.phi_f_1, x, P.phi_s_2, P.phi_r_2);
        %         fun = @(x) P.phi0*ones(size(x));
    case 'rho'
        fun = @(x) sigmoid_dec(P.rho_f_0, P.rho_f_1, x, P.rho_s_2, P.rho_r_2);
        %         fun = @(x)  P.rho0*ones(size(x));
    case 'psi'
        fun = @(x) sigmoid_dec(P.psi_f_0, P.psi_f_1, x, P.psi_s_2, P.psi_r_2);
         %         fun = @(x) P.psi0*ones(size(x));
    otherwise
        error('not defined probability parameter')
end

end

% f_0  value min
% f_1  value max
% t_2  shift
% s_2  sigmoid steepness, smaller is steeper


