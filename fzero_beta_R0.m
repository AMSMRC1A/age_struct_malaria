%% find the value of betaM that gives a target R0 value (used for generating cuts in the bifurcation plots)
clear all
close all
clc
format long
global P
%% numerical config
tfinal = 100*365; age_max = 80*365; P.age_max = age_max;
dt = 100; % time/age step size in days, default = 50; could go dt = 200 (still robust)
da = dt; t = (0:dt:tfinal)'; nt = length(t); a = (0:da:age_max)'; na = length(a);

P.a = a; P.na = na; P.nt = nt; P.dt = dt; P.da = da; P.t = t;

Malaria_parameters_baseline; % model parameters - rates are in 1/day
lP = 'betaM';  % bifurcating parameters
immunity_feedback = -1; % -1 = low immunity fixed, 1 = dynamic, 0 = high fixed
if immunity_feedback == -1 % betaM = 0.008
    P.phi_f_0 = 0.230971153687268; % value at zero
    P.phi_f_1 = 0.230971153687268; % value at L (function saturates to this value)
    P.rho_f_0 = 0.907631179204690; % value at zero
    P.rho_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
    P.psi_f_0 = 0.907631179204690 ; % value at zero
    P.psi_f_1 = 0.907631179204690 ; % value at L (function saturates to this value)
    % 0.138 gives R0=5
    % 0.0885 gives R0=4
elseif immunity_feedback == 0 % betaM = 0.25
    % average populational sigmoids f0 = f1 = average
    P.phi_f_0 = 0.915792480087329; % value at zero
    P.phi_f_1 = 0.915792480087329; % value at L (function saturates to this value)
    P.rho_f_0 = 0.114825053290306; % value at zero
    P.rho_f_1 = 0.114825053290306; % value at L (function saturates to this value)
    P.psi_f_0 = 0.114825053290306; % value at zero
    P.psi_f_1 = 0.114825053290306; % value at L (function saturates to this value)
    % 0.538872779193500 gives R0=5
    % 0.344878578683840 gives R0=4  
elseif immunity_feedback == 1
    % 0.539 gives R0=5
    % 0.0833  gives R0=4
end

% R0_target = 4;
R0_target = 1;
beta_0 = [0 1];
func = @(beta) R0_beta(beta)-R0_target; 
beta_zero = fzero(func,beta_0)
round(beta_zero,3,'significant')

function y = R0_beta(beta)
global P
P.betaM = beta;
y = R0_cal;
end