function [bH,bM] = biting_rate(NH,NM)
global P
% happy mosquitoes:
% bM = P.bm;
% bH = P.bm*NM/NH;

% compromised model
b_tot = P.bm*P.bh*NH*NM/(P.bm*NM+P.bh*NH);
bH = b_tot/NH; % bites per human
bM = b_tot/NM; % bites per mosquito

end