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
        % Cac
        fun =  @(alpha) integral(@(a) exp(a/P.dac).*v(a),0,alpha);
        Cac = @(alpha) P.cV*exp(-alpha/P.dac).*(v(0)+fun(alpha));
        % Cm
        fun = @(alpha) exp(-alpha/P.dac).*gH(alpha).*exp(-muH_int(alpha));
        Cm01 = P.cV*P.m*v(0)*integral(@(alpha) fun(alpha),0,Inf);
        amax = @(alpha) alpha;
        Cm02 = P.cV*P.m*integral2(@(alpha,a) fun(alpha).*exp(a./P.dac).*v(a), 0, inf, 0, amax);
        Cm0 = Cm01+Cm02;
        Cm = @(alpha) Cm0*exp(-alpha./P.dm);
        % CH
        CH = @(alpha) P.c1*Cac(alpha)+P.c2*Cm(alpha);     
    otherwise
        error('undefined steady state label...')
end
end