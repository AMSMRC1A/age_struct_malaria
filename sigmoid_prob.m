%% sigmoidal function for converting immunity to probabilties
function f = sigmoid_prob(x, lprob)

sigmoid_fun = sigmoid_prob_fun(lprob); % return a function handle
f = sigmoid_fun(x);

end