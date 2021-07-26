function [bH,bM] = biting_rate(NH,NM)
global P
% NH(1,np) NM(1,np) - row vectors
% b_tot bH bM dimension 1*np - row vectors

% compromise model from Chitniss et al. (2006)
b_tot = P.bm*P.bh*NH.*NM./(P.bm*NM+P.bh*NH);
bH = b_tot./NH; % bites per human
bM = b_tot./NM; % bites per mosquito

end

% happy mosquitoes:
% bM = P.bm;
% bH = P.bm*NM/NH;