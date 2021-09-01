%% sigmoidal function for converting immunity to probabilties
function f = sigmoid_prob(x, lprob)
global P

switch lprob
    case 'phi'
        cmin = 0; cmax = 1; k = 2; c = P.halfmaximal; % c = 40 (if in years, i.e. da = 365, Filipe)
    case 'rho'
        cmin = 0; cmax = 1; k = 2; c = P.halfmaximal;
    case 'psi'
        cmin = 0; cmax = 1; k = 2; c = P.halfmaximal;
    otherwise 
        error('not defined probability parameter')
end
    f = 0.5*ones(size(x));
    %f = cmax*((1-cmin)*c^k./(c^k+x.^k)+cmin);
    %if lprob == 'phi'
    %    f = 1 - f;
    %end
end

% input x should be the immunity level (C_s)
% minimum value: cmin
% maximum value: cmax
% constant at half maximal: c
% Hill coefficient (steepness): k