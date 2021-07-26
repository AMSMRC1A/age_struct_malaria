function [SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM)
% evolve mosquito quantities, keep at the steady state 
% NH(1,np) NM(1,np) - row vectors
% SM(1,np)
% lamM(1,np)  bM(1,np)
% DH NH_star(na,np)
global P

na = P.na;
np = P.np;

DH = reshape(DH,na,np);
AH = reshape(AH,na,np);

NH_star = NH*P.pmat'; % dimension 1*np
DH_star = DH*P.pmat'; % dimension na*np
AH_star = AH*P.pmat'; % dimension na*np

lamM = FOI_M(DH_star,AH_star,NH_star,NM);

SM = P.gM./(lamM+P.muM);
IM = (P.gM/P.muM)*(P.sigma/(P.sigma+P.muM))*(lamM./(lamM+P.muM));
EM = (P.muM/P.sigma)*IM;

end