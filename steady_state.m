function [S,E,D,A,Cac,Cm,Ctot] = steady_state(lstate)
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
        Ctot = @(alpha) P.c1*Cac(alpha)+P.c2*Cm(alpha);
    case 'EE'
        % ***only tested for no immunity feedback***
        % use numerical simulation for an initial guess
        dt = 20; tfinal= 10*365; % run for a few years to get closer to EE
        t = (0:dt:tfinal)'; nt = length(t);
        tic
        [SH, EH, DH, AH, ~, ~, ~, ~, Cac, ~] = age_structured_Malaria(P.na,P.da,nt);
        toc                
        %% run solver on a coarser grid to speed up
        da_fine = P.da; da_coarse = 100; P.da = da_coarse; 
        a_fine = P.a; a_coarse = (0:da_coarse:P.age_max)'; P.a = a_coarse;
        na_fine = P.na; na_coarse = length(a_coarse); P.na = na_coarse;
        Malaria_parameters_transform;
        SH0 = interp1(a_fine,SH(:,end),a_coarse)./P.PH_stable; 
        EH0 = interp1(a_fine,EH(:,end),a_coarse)./P.PH_stable;
        DH0 = interp1(a_fine,DH(:,end),a_coarse)./P.PH_stable;
        AH0 = interp1(a_fine,AH(:,end),a_coarse)./P.PH_stable;
        Cac0 = interp1(a_fine,Cac(:,end),a_coarse)./P.PH_stable;
        x0 = [SH0;EH0;DH0;AH0;Cac0];
        options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
        F_prop = @(x) human_model_der_prop(x);
        tic
        [xsol,err,~,~,~] = fsolve(F_prop,x0,options);
        toc
        if max(max(abs(err)))>10^-5
            disp('not converged')
            keyboard
        end       
        x_EE = reshape(xsol,[P.na,5]);
        S_coarse = x_EE(:,1);
        E_coarse = x_EE(:,2);
        D_coarse = x_EE(:,3);
        A_coarse = x_EE(:,4);
        Cac_coarse = x_EE(:,5);      
        %% return the values on fine grid with interpolation
        % recover numerical config
        P.da = da_fine; P.a = a_fine; P.na = na_fine;
        Malaria_parameters_transform;
        S = interp1(a_coarse,S_coarse,a_fine).*P.PH_stable; 
        E = interp1(a_coarse,E_coarse,a_fine).*P.PH_stable; 
        D = interp1(a_coarse,D_coarse,a_fine).*P.PH_stable; 
        A = interp1(a_coarse,A_coarse,a_fine).*P.PH_stable; 
        Cac = interp1(a_coarse,Cac_coarse,a_fine).*P.PH_stable; 
        Cm0 = P.m*trapz(P.gH.*Cac)*P.da/P.PH_stable(1);
        Cm = Cm0.*exp(-P.a./P.dm).*P.PH_stable;
        Ctot = P.c1*Cac+P.c2*Cm;
    otherwise
        error('undefined steady state label...')
end
end