function Q_val = QOI_value(lQ)
global lP
global P flag_disp

FileName = ['Results/SA/',lQ,'_',lP,'_',num2str(P.(lP),5),'.mat'];

if exist(FileName,'file') % Q_val already calculated before
    if flag_disp; disp('load Q_val results...'); end
    load(FileName,'Q_val')
    return
end

switch lQ
    case 'R0'    
        Q_val  = R0_cal();
    otherwise
        keyboard
end

save(FileName,'Q_val')

end
