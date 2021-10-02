function [S,E,D,A,Cac,Cm,CH] = steady_state(lstate)
% return function handles for 'DFE', numerical for 'EE'
% all variables are pop size
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
        % ***only tested for no immunity feedback***
        % use rough numerical simulation for an initial guess
        dt_new = 20; da_new = dt_new; 
        tfinal_new = 10*365; % run for 10 years to get closer to EE
        t_new = (0:dt_new:tfinal_new)'; nt_new = length(t_new);
        a_new = (0:da_new:P.age_max)'; na_new = length(a_new);
        a_old = P.a; na_old = P.na; nt_old = P.nt; dt_old = P.dt; da_old = P.da; t_old = P.t; % save
        P.a = a_new; P.na = na_new; P.nt = nt_new; P.dt = dt_new; P.da = da_new; P.t = t_new;
        tic
        [SH, EH, DH, AH, ~, ~, ~, ~, ~, ~] = age_structured_Malaria();
        toc
        P.a = a_old; P.na = na_old; P.nt = nt_old; P.dt = dt_old; P.da = da_old; P.t = t_old; % recover
        x0 = [SH(:,end)./P.PH_stable;EH(:,end)./P.PH_stable;DH(:,end)./P.PH_stable;AH(:,end)./P.PH_stable];
        options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
        F_prop = @(x) human_model_der_prop(x);
        tic
        [xsol,err,~,~,~] = fsolve(F_prop,x0,options);
        toc
        if max(max(abs(err)))>10^-5
            disp('not converged')
            keyboard
        end
        x_EE = reshape(xsol,[P.na,4]);
        S = x_EE(:,1).*P.PH_stable;
        E = x_EE(:,2).*P.PH_stable;
        D = x_EE(:,3).*P.PH_stable;
        A = x_EE(:,4).*P.PH_stable;
        Cac = NaN(size(S)); % not implemented
        Cm = NaN(size(S));  % not implemented
        CH = P.c1*Cac+P.c2*Cm; % not implemented       
    otherwise
        error('undefined steady state label...')
end
end