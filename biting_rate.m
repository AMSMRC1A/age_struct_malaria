function [bH,bM] = biting_rate(NH,NM)
global P
% happy mosquitoes:
bM = P.bm;
bH = P.bm*NM/NH;

end