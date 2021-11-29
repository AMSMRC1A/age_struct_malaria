function [S,E,D,A,Cac,Cm,Ctot] = steady_state(lstate,lreturn)
% lreturn = 'handle' return function handles for 'DFE'
% lreturn = 'numerical' return numerical values
% all variables are pop size
global P

gH = P.gH_fun; % feritlity function handle
muH_int = P.muH_int_fun; % death function handle M
switch lstate
    case 'DFE'
        if strcmp(lreturn,'handle')
            S = @(alpha) 1*ones(size(alpha)).*P.PH_stable_fun(alpha);
            E = @(alpha) 0*ones(size(alpha)).*P.PH_stable_fun(alpha);
            D = @(alpha) 0*ones(size(alpha)).*P.PH_stable_fun(alpha);
            A = @(alpha) 0*ones(size(alpha)).*P.PH_stable_fun(alpha);
            % Cac exact expression obtained based on simple v(alpha)
            Cac_prop = @(alpha) P.cV*exp(-alpha/P.dac).*(P.vb0+P.vb0*P.dac*(exp(alpha/P.dac)-1));
            % Cm
            Cm0 = P.m*integral(@(alpha) gH(alpha).*exp(-muH_int(alpha)).*Cac_prop(alpha), 0, P.age_max);
            Cm_prop = @(alpha) Cm0.*exp(-alpha./P.dm);
            % CH
            Cac = @(alpha) Cac_prop(alpha).*P.PH_stable_fun(alpha);
            Cm = @(alpha) Cm_prop(alpha).*P.PH_stable_fun(alpha);
            Ctot = @(alpha) P.c1*Cac(alpha)+P.c2*Cm(alpha);
        elseif strcmp(lreturn,'numerical')
            a = P.a;
            S = 1*ones(size(a)).*P.PH_stable;
            E = 0*ones(size(a)).*P.PH_stable;
            D = 0*ones(size(a)).*P.PH_stable;
            A = 0*ones(size(a)).*P.PH_stable;
            % Cac exact expression obtained based on simple v(alpha)
            Cac_prop = P.cV*exp(-a./P.dac).*(P.vb0+P.vb0*P.dac*(exp(a./P.dac)-1));
            Cac = Cac_prop.*P.PH_stable;
            Cm0 = P.m*trapz(P.gH.*exp(-P.muH_int).*Cac_prop)*P.da;
            Cm = Cm0.*exp(-a./P.dm).*P.PH_stable;
            Ctot = P.c1*Cac+P.c2*Cm;
        end
    case 'EE'
        if strcmp(lreturn,'numerical')
            dt = P.dt; tfinal= 50*365;  % run for a long time; numerical EE
            da = dt; a = (0:da:P.age_max)'; na = length(a);
            [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
            [SH, EH, DH, AH, ~, ~, ~, Cmt, Cact, Ctott] = age_structured_Malaria(da,na,tfinal, SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
            S = SH(:,end); E = EH(:,end); D = DH(:,end); A = AH(:,end); 
            Cac = Cact(:,end); Cm = Cmt(:,end); Ctot = Ctott(:,end);
        elseif strcmp(lreturn,'fsolve')
            % use numerical simulation for an initial guess
            dt = P.dt; tfinal= 15*365;  % run for a few years to get closer to EE
            da = dt; a = (0:da:P.age_max)'; na = length(a);
            [SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0] = age_structured_Malaria_IC('init');
            [SH, EH, DH, AH, ~, ~, ~, ~, Cac, ~] = age_structured_Malaria(da,na,tfinal, SH0, EH0, DH0, AH0, SM0, EM0, IM0, Cm0, Cac0, Ctot0);
            %% run solver on a coarser grid to speed up
            da_fine = P.da; da_coarse = 80; P.da = da_coarse;
            a_fine = P.a; a_coarse = (0:da_coarse:P.age_max)'; P.a = a_coarse;
            na_fine = P.na; na_coarse = length(a_coarse); P.na = na_coarse;
            Malaria_parameters_transform;
            SH0 = interp1(a_fine,SH(:,end),a_coarse);
            EH0 = interp1(a_fine,EH(:,end),a_coarse);
            DH0 = interp1(a_fine,DH(:,end),a_coarse);
            AH0 = interp1(a_fine,AH(:,end),a_coarse);
            Cac0 = interp1(a_fine,Cac(:,end),a_coarse);
            PH0 = SH0 + EH0 + AH0 + DH0;
            x0 = [SH0./PH0; EH0./PH0; DH0./PH0; AH0./PH0; Cac0./PH0];
            options = optimoptions('fsolve','Display','none','OptimalityTolerance', 1e-25);
            F_prop = @(x) human_model_der_prop(x);
            [xsol,err,~,~,~] = fsolve(F_prop,x0,options);
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
            % update progression probability based on immunity Ctot
            P.phi = sigmoid_prob(Ctot./P.PH_stable, 'phi'); % prob. of DH -> RH
            P.rho = sigmoid_prob(Ctot./P.PH_stable, 'rho'); % prob. of severely infected, EH -> DH
            P.psi = sigmoid_prob(Ctot./P.PH_stable, 'psi'); % prob. AH -> DH
        end
    otherwise
        error('undefined steady state label...')
end
end