function Qp_val = QOIp_value(lQ,lP)
      
global P flag_disp 

factor = 0.001;

if flag_disp; disp(['dp=',num2str(factor)]); end

switch lP
    case {'v0','bh','bm','betaM','betaD','betaA','muM','MHr','sigma','phi_f_0', 'phi_f_1', 'phi_t_2', 'phi_s_2', 'rho_f_0', 'rho_f_1', 'rho_t_2', 'rho_s_2', 'psi_f_0', 'psi_f_1', 'psi_t_2', 'psi_s_2'} % use local numerical derivatives
        P_base = P.(lP);  
        delta_P = abs(P_base * factor);
        P.(lP) = P_base-delta_P; % p_lower 
        Malaria_parameters_transform;
        if flag_disp; disp('P-dp run ---'); end
        Q_lower = QOI_value(lQ);
        P.(lP) = P_base+delta_P; % p_upper
        Malaria_parameters_transform;    
        if flag_disp; disp('P+dp run ---'); end
        Q_upper = QOI_value(lQ);
        Qp_val = (Q_upper-Q_lower)/(2*delta_P); 
    otherwise
        keyboard
end

end 