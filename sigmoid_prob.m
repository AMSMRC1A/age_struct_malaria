%% sigmoidal function for converting immunity to probabilties

function f = sigmoid_prob(x, cmin, cmax, c, k)
    f = cmax*((1-cmin)*c^k./(c^k+x.^k)+cmin);
end

% input x should be the immunity level (C_s)
% minimum value: cmin
% maximum value: cmax
% constant at half maximal: c
% Hill coefficient (steepness): k