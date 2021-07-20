function lamM = FOI_M(bM,DH,AH,NH)
global P
% DH and AH are vectors at time t=n;
% bM is calculated in biting_rate.m

da = P.da;
lamM = bM*trapz(P.betaD*DH+P.betaA*AH)*da/NH; 

end