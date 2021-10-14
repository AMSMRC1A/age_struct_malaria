%% sigmoidal function for converting immunity to probabilties
% return function handle
function fun = sigmoid_prob_fun(lprob)
global P

sigmoid_inc = @(f_0, f_1, x, t_2, s_2, L) f_0 + (f_1-f_0)./(1 + exp(-(x/L-t_2)/s_2));
sigmoid_dec = @(f_0, f_1, x, t_2, s_2, L) f_1 + (f_0-f_1)./(1 + exp(-(x/L-t_2)/s_2));

% NB set f1 = f0 to get a constant function (equal to f0 everywhere)

switch lprob
    case 'phi'
        fun = @(x) sigmoid_inc(P.phi_f_0, P.phi_f_1, x, P.phi_t_2, P.phi_s_2, P.L);
    case 'rho'
        fun = @(x) sigmoid_dec(P.rho_f_0, P.rho_f_1, x, P.rho_t_2, P.rho_s_2, P.L);
    case 'psi'
        fun = @(x) sigmoid_dec(P.psi_f_0, P.psi_f_1, x, P.psi_t_2, P.psi_s_2, P.L);
    otherwise
        error('not defined probability parameter')
end

end

% L effective range is [0,L]
% f_0  value at zero
% f_1  value at L (function saturates to this value)
% t_2  threshold value (as a fraction of L)
% s_2  sigmoid steepness, smaller is steeper