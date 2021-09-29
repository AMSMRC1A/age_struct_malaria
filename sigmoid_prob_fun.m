%% sigmoidal function for converting immunity to probabilties
% return function handle
function fun = sigmoid_prob_fun(lprob)
global P

switch lprob
    case 'phi'
        cmin = 0; cmax = 1; k = 2; c = 40; % c = 40 (if in years, i.e. da = 365, Filipe)
        fun = @(x) P.phi0*ones(size(x));
    case 'rho'
        cmin = 0; cmax = 1; k = 2; c = 40;
        fun = @(x) P.rho0*ones(size(x));
    case 'psi'
        cmin = 0; cmax = 1; k = 2; c = 40;
        fun = @(x) P.psi0*ones(size(x));
    otherwise
        error('not defined probability parameter')
end
% fun = @(x) cmax*((1-cmin)*c^k./(c^k+x.^k)+cmin);
end

% input x should be the immunity level (C_s)
% minimum value: cmin
% maximum value: cmax
% constant at half maximal: c
% Hill coefficient (steepness): k