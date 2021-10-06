function Q_val = QOI_value(lQ)
global lP
global P flag_disp
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

da = P.da;
a = P.a;
use_EE_data = 0;

if strcmp(lQ(1:2),'EE')
    FileName_EE = ['Results/SA/EE_',lP,'_',num2str(P.(lP),5),'.mat'];
    if exist(FileName_EE,'file') && use_EE_data==1 % EE calculated before
%         disp('load EE results...');
        EE = load(FileName_EE,'SH','EH','DH','AH','P');
        AA = EE.P; % load
        BB = P; % current
        fields = {'v_fun','muH_int_fun','gH_fun'};
        AA = rmfield(AA,fields); BB = rmfield(BB,fields);
        [~,d1,d2] = comp_struct(AA,BB,0,0,10^-16);
        if ~isempty(d1) || ~isempty(d2)
            disp('baseline values changed, need to recal EE: remove old EE files to other places')
            keyboard
        end
        SH = EE.SH; EH = EE.EH; DH = EE.DH; AH = EE.AH;
    else
        if R0_cal()<1
            disp('EE DNE')
            keyboard
        end
        [SH,EH,DH,AH,Cac,Cm,Ctot] = steady_state('EE');
        %%
        figure_setups; hold on;
        plot(a/365,SH,'-','Color',colour_mat1);
        plot(a/365,EH,'-','Color',colour_mat3);
        plot(a/365,DH,'-','Color',colour_mat2);
        plot(a/365,AH,'-','Color',colour_mat7);
        plot(a/365,SH+EH+DH+AH,'-k');
        legend('SH (solver)','EH (solver)','DH (solver)', 'AH (solver)', 'PH (solver)');
        title('Final Age Dist.');
        xlabel('age');
        grid on
        axis([0 P.age_max/365 0 max(SH+EH+DH+AH)]);
        
        figure_setups; hold on;
        plot(a/365,Cac,'-.r');
        plot(a/365,Cm,'-.b');
        plot(a/365,Ctot,'-.k');
        xlabel('age (years)')
        legend('Acquired','Maternal','Total','Location','SouthEast');
        title('Immun dist.');
        axis([0 P.age_max/365 0 max(Ctot)*1.1]);
        grid on
        keyboard
%         save(FileName_EE,'SH','EH','DH','AH','P')
    end
end

switch lQ
    case 'R0'    
        Q_val  = R0_cal();
    case 'RHM'
        [~,Q_val,~]  = R0_cal();
    case 'RMH'
        [~,~,Q_val]  = R0_cal();
    case 'EE-infected'
        Q_val = 1-da*trapz(SH);       
    case 'EE-DA'
        Q_val = 1-da*trapz(SH)-da*trapz(EH);
    case 'EE-D-frac'
        Q_val = trapz(DH)/trapz(DH+AH);
    case 'EE-EDA'
        NH = trapz(SH+EH+DH+AH)*da; 
        Q_val = [trapz(EH)*da/NH; trapz(DH)*da/NH; trapz(AH)*da/NH];
    case 'EE-EIR'
        NH = trapz(SH+EH+DH+AH)*da; 
        NM = P.gM/P.muM;
        [bH,bM] = biting_rate(NH,NM);
        Lambda_M = bM*trapz(P.betaD*DH + P.betaA*AH)*da;
        IM_frac_EE = P.sigma/(P.sigma+P.muM)*(Lambda_M/(Lambda_M + P.muM));   
        Q_val = bH*IM_frac_EE*365; % annual EIR           
    otherwise
        keyboard
end



end
