function Q_val = QOI_value(lQ)
global lP
global P flag_disp
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

age_max = P.age_max;
a = P.a;
t = P.t;
da = P.da;
 
FileName = ['Results/SA/',lQ,'_',lP,'_',num2str(P.(lP),5),'.mat'];

if exist(FileName,'file') % Q_val already calculated before
    if flag_disp; disp('load Q_val results...'); end
    load(FileName,'Q_val')
    return
end

switch lQ
    case 'R0'    
        Q_val  = R0_cal();
    case 'RHM'
        [~,Q_val,~]  = R0_cal();
    case 'RMH'
        [~,~,Q_val]  = R0_cal();
    case 'EIR_EE'
        %% EIR
        [SH, EH, DH, AH, SM, EM, IM, ~, ~, ~] = age_structured_Malaria;
        PH_final = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end); % total human at age a, t = n
        NH = trapz(PH_final)*da;
        NM = SM(end)+EM(end)+IM(end);
        [bh,~] = biting_rate(NH,NM);
        Q_val = bh.*IM(end)./NM*365; % annual EIR
    case 'EE_numerical'
        % use numerical simulation
        [SH, EH, DH, AH, SM, EM, IM, ~, ~, ~] = age_structured_Malaria;
        PH = SH+EH+DH+AH; % total human at age a, t = n
        F_prop = @(x) human_model_der_fun(x);
        F_number = @(x) human_model_der_fun2(x);
        x_EE = [SH(:,end);EH(:,end);DH(:,end);AH(:,end)];
        if norm(F_number(x_EE)) > 10^-5
            disp('EE not achieved')
            keyboard
        end
        max(max(abs(F_number(x_EE))))     
        x_EE_prop = [SH(:,end)./PH(:,end);EH(:,end)./PH(:,end);DH(:,end)./PH(:,end);AH(:,end)./PH(:,end)];
        max(max(abs(F_prop(x_EE_prop))))  
        keyboard
        %% plot
        figure_setups; hold on;
        plot(a,SH(:,end),'-','Color',colour_mat1);
        plot(a,EH(:,end),'-','Color',colour_mat3);
        plot(a,DH(:,end),'-','Color',colour_mat2);
        plot(a,AH(:,end),'-','Color',colour_mat7);
        plot(a,PH(:,end),'-k');
        legend('SH (numercial)','EH (numercial)','AH (numercial)', 'DH (numercial)','PH (numercial)');
        title('Final Age Dist.');
        xlabel('age');
        axis_years(gca,age_max); % change to x-axis to years if needed
        axis([0 age_max 0 max(SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end))]);
        grid on; grid minor
        figure_setups; hold on;
        plot(a,SH(:,end)./PH(:,end),'-','Color',colour_mat1);
        plot(a,EH(:,end)./PH(:,end),'-','Color',colour_mat3);
        plot(a,DH(:,end)./PH(:,end),'-','Color',colour_mat2);
        plot(a,AH(:,end)./PH(:,end),'-','Color',colour_mat7);
        plot(a,PH(:,end)./PH(:,end),'-k');
        legend('SH prop (numercial)','EH prop (numercial)','DH prop (numercial)', 'AH prop (numercial)','PH (numercial)');
    case 'EE_fsolve'
        error_tolerance = 1e-25;
        options = optimoptions('fsolve','Display','none','OptimalityTolerance', error_tolerance);
        F_prop = @(x) human_model_der_fun(x);
        x0 = [1; 0.4*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1); 0; 0.2*ones(length(a)-1,1)]; % initial guess for the EE
        [xsol,error,~,~,jacobian] = fsolve(F_prop,x0,options);

        
%         F_number = @(x) human_model_der_fun2(x);
%         x_EE = reshape(xsol,[P.na,4]);
%         figure_setups; hold on;
%         plot(a,x_EE(:,1),'-','Color',colour_mat1); 
%         plot(a,x_EE(:,2),'-','Color',colour_mat3);
%         plot(a,x_EE(:,3),'-','Color',colour_mat2);
%         plot(a,x_EE(:,4),'-','Color',colour_mat7);
%         plot(a,sum(x_EE,2),'-k');
%         axis([0 age_max 0 1]);
%         legend('SH prop (solver)','EH prop (solver)','DH prop (solver)', 'AH prop (solver)', 'PH prop (solver)');
%         title('Final Age Dist.');
%         xlabel('age');
%         axis_years(gca,P.age_max); % change to x-axis to years if needed
%         grid on
%         x_EE1(:,1) = x_EE(:,1).*P.PH_stable;
%         x_EE1(:,2) = x_EE(:,2).*P.PH_stable;
%         x_EE1(:,3) = x_EE(:,3).*P.PH_stable;
%         x_EE1(:,4) = x_EE(:,4).*P.PH_stable;
%         figure_setups; hold on;
%         plot(a,x_EE1(:,1),'-','Color',colour_mat1); 
%         plot(a,x_EE1(:,2),'-','Color',colour_mat3);
%         plot(a,x_EE1(:,3),'-','Color',colour_mat2);
%         plot(a,x_EE1(:,4),'-','Color',colour_mat7);
%         plot(a,sum(x_EE1,2),'-k');        
%         axis([0 age_max 0 max(sum(x_EE1,2))]);
%         legend('SH prop (solver)','EH prop (solver)','DH prop (solver)', 'AH prop (solver)', 'PH prop (solver)');
%         title('Final Age Dist.');
%         xlabel('age');
%         axis_years(gca,P.age_max); % change to x-axis to years if needed
%         grid on
%         keyboard
        
%         x0 = reshape(x_EE1,P.na*4,1);
%         [xsol,error,~,~,~] = fsolve(F_number,x0,options);
%         max(max(abs(error)))
%         x_EE2 = reshape(xsol,[P.na,4]);
%         figure_setups; hold on;
%         plot(a,x_EE2(:,1),'-','Color',colour_mat1); 
%         plot(a,x_EE2(:,2),'-','Color',colour_mat3);
%         plot(a,x_EE2(:,3),'-','Color',colour_mat2);
%         plot(a,x_EE2(:,4),'-','Color',colour_mat7);
%         plot(a,sum(x_EE2,2),'-k');
%         axis([0 age_max 0 max(sum(x_EE2,2))]);
%         legend('SH prop (solver)','EH prop (solver)','DH prop (solver)', 'AH prop (solver)', 'PH prop (solver)');
%         title('Final Age Dist.');
%         xlabel('age');
%         axis_years(gca,P.age_max); % change to x-axis to years if needed
%         grid on
%         keyboard
%         x_EE(:,1) = x_EE(:,1)./P.PH_stable;
%         x_EE(:,2) = x_EE(:,2)./P.PH_stable;
%         x_EE(:,3) = x_EE(:,3)./P.PH_stable;
%         x_EE(:,4) = x_EE(:,4)./P.PH_stable;
%         rr = F_prop(reshape(x_EE,P.na*4,1))
%         max(max(abs(rr)))
%         keyboard
    case 'stability'
        
    otherwise
        keyboard
end

% save(FileName,'Q_val')

end
