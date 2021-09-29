function Q_val = QOI_value(lQ)
global lP
global P flag_disp
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

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
    case 'EE_infected'
        [SH,~,~,~,~,~,~] = steady_state('EE');
        Q_val = 1-da*trapz(SH(:,1).*P.PH_stable);
    case 'EE_EDA'
        [SH,EH,DH,AH,~,~,~] = steady_state('EE');
        NH = trapz(SH+EH+DH+AH)*da; 
        Q_val = [trapz(EH)*da/NH; trapz(DH)*da/NH; trapz(AH)*da/NH];
        save(FileName,'Q_val')
    case 'EE_EIR'
        [SH,EH,DH,AH,~,~,~] = steady_state('EE');
        NH = trapz(SH+EH+DH+AH)*da; 
        NM = P.gM/P.muM;
        [bH,bM] = biting_rate(NH,NM);
        Lambda_M = bM*trapz(P.betaD*DH + P.betaA*AH)*da;
        IM_frac_EE = P.sigma/(P.sigma+P.muM)*(Lambda_M/(Lambda_M + P.muM));   
        Q_val = bH*IM_frac_EE*365; % annual EIR           
    otherwise
        keyboard
end


% save(FileName,'Q_val')

end
