function [S,E,D,A,Cac,Cm,CH] = steady_state(lstate)
% return function handles for steady states
global P

gH = P.gH_fun; % feritlity function handle
v = P.v_fun; % vaccination function handle
muH_int = P.muH_int_fun; % death function handle M

switch lstate
    case 'DFE'
        S = @(alpha) 1;
        E = @(alpha) 0;
        D = @(alpha) 0;
        A = @(alpha) 0;        
        % Cac  exact expression obtained based on simple v(alpha)
        Cac = @(alpha) P.cV*exp(-alpha/P.dac).*(P.v0+P.v0*P.dac*(exp(alpha/P.dac)-1));
        % Cm
        Cm0 = P.m*integral(@(alpha) gH(alpha).*exp(-muH_int(alpha)).*Cac(alpha), 0, P.age_max);        
        Cm = @(alpha) Cm0.*exp(-alpha./P.dm);        
        % CH
        CH = @(alpha) P.c1*Cac(alpha)+P.c2*Cm(alpha);     
    otherwise
        error('undefined steady state label...')
end
end