%% sigmoidal function for converting immunity to probabilties
function f = sigmoid_prob(x, lprob)
% x(na,np) f(na,np)
global P

x = reshape(x,P.na,P.np);

switch lprob
    case 'phi'
        cmin = 0; cmax = 1; k = 2; c = 40; % c = 40 (if in years, i.e. da = 365, Filipe)
    case 'rho'
        cmin = 0; cmax = 1; k = 2; c = 40;
    case 'psi'
        cmin = 0; cmax = 1; k = 2; c = 40;
    otherwise 
        error('not defined probability parameter')
end
    %f = cmax*((1-cmin)*c^k./(c^k+x.^k)+cmin);
    f = 0.5*ones(size(x));
end

% input x should be the immunity level (C_s)
% minimum value: cmin
% maximum value: cmax
% constant at half maximal: c
% Hill coefficient (steepness): k