function y = f(x)
% function on FOI in the Cac immunity equation
global P
   y = x/(P.uc*x+1); % saturated function used in Griffin et al (2015)
end