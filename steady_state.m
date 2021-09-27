function [S,E,D,A,Cac,Cm,CH] = steady_state(lstate)
% return function handles for 'DFE', numerical for 'EE'
global P

gH = P.gH_fun; % feritlity function handle
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
    case 'EE'
        % use fsolve for RHS
        error_tolerance = 1e-20;
        options = optimoptions('fsolve','Display','none','OptimalityTolerance', error_tolerance);
        F = @(x) human_model_der_fun(x);
        x0 = [1; 0.4*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1)]; % initial guess for the EE
        [xsol,~,~,~,~] = fsolve(F,x0,options);
        x_EE = reshape(xsol,[P.na,4]);
        S = x_EE(:,1);
        E = x_EE(:,2);
        D = x_EE(:,3);
        A = x_EE(:,4);
        Cac = NaN(size(S));
        Cm = NaN(size(S));
        CH = P.c1*Cac+P.c2*Cm;
        % use numerical steady state
        
    otherwise
        error('undefined steady state label...')
end
end