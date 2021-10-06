%% sigmoidal function for converting immunity to probabilties
% return function handle
function fun = sigmoid_prob_fun(lprob)
global P

sigmoid = @(f_0, f_1, x, t_2, s_2, L) f_0 + (f_1-f_0)./(1 + exp(-(x/L-t_2)/s_2));

% NB set f1 = f0 to get a constant function (equal to f0 everywhere)

switch lprob
    case 'phi'
        %cmin = 0; cmax = 1; k = 2; c = 40; % c = 40 (if in years, i.e. da = 365, Filipe)
        fun = @(x) sigmoid(P.phi_f_0, P.phi_f_1, x, P.phi_t_2, P.phi_s_2, P.L);
        %P.phi0*ones(size(x));
    case 'rho'
        %cmin = 0; cmax = 1; k = 2; c = 40;
        fun = @(x) sigmoid(P.rho_f_0, P.rho_f_1, x, P.rho_t_2, P.rho_s_2, P.L);
        %P.rho0*ones(size(x));
    case 'psi'
        %cmin = 0; cmax = 1; k = 2; c = 40;
        fun = @(x) sigmoid(P.psi_f_0, P.psi_f_1, x, P.psi_t_2, P.psi_s_2, P.L);
        %P.psi0*ones(size(x));
    otherwise
        error('not defined probability parameter')
end
end

% input x should be the immunity level (C_s)
% minimum value: cmin
% maximum value: cmax
% constant at half maximal: c
% Hill coefficient (steepness): k