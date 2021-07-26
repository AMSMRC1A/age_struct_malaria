function lamM = FOI_M(DH,AH,NH,NM)
global P
% DH and AH are vectors at time t=n;
% bM is calculated in biting_rate.m

% NH(1,np) AH(na,np)
% trapz integrates over each column and returns a row vector
% lamM(1,np) bM(1,np)

da = P.da;
[~,bM] = biting_rate(NH,NM); 
lamM = bM.*trapz(P.betaD*DH+P.betaA*AH)*da./NH; 

end