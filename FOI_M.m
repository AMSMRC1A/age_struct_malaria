function lamM = FOI_M(bM,DH,AH,NH)
global P
da = P.da;

% DH and AH are vectors at time t=n;
lamM = bM*trapz(P.betaD*DH+P.betaA*AH)*da/NH; 

end