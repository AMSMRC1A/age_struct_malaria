%% sigmoidal function for converting immunity to probabilties
function f = sigmoid_prob(x, lprob)
switch lprob
    case 'phi'
        cmin = 0; cmax = 1; k = 2; c = 0.0001;
    case 'rho'
        cmin = 0; cmax = 1; k = 2; c = 0.0001;
    case 'psi'
        cmin = 0; cmax = 1; k = 2; c = 0.0001;
    otherwise 
        error('not defined probability parameter')
end
    f = cmax*((1-cmin)*c^k./(c^k+x.^k)+cmin);
end

% input x should be the immunity level (C_s)
% minimum value: cmin
% maximum value: cmax
% constant at half maximal: c
% Hill coefficient (steepness): k