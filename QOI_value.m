function Q_val = QOI_value(lQ)
global lP
global P flag_disp
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

da = P.da;

if strcmp(lQ(1:2),'EE')
    FileName_EE = ['Results/SA/EE_',lP,'_',num2str(P.(lP),5),'.mat'];
    if exist(FileName_EE,'file') % EE calculated before
        disp('load EE results...');
        EE = load(FileName_EE,'SH','EH','DH','AH');
        SH = EE.SH; EH = EE.EH; DH = EE.DH; AH = EE.AH;
    else
        [SH,EH,DH,AH,~,~,~] = steady_state('EE');
        save(FileName_EE,'SH','EH','DH','AH')
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
