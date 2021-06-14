function lamM = FOI_M(bM,DH,AH,NH)
global P
% DH and AH are vectors at time t=n;

da = P.da;
lamM = bM*trapz(P.betaD*DH+P.betaA*AH)*da/NH; 

end