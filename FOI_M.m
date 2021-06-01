function lamM = FOI_M(bM,DH,AH,NH)
global P
global a 
% DH and AH are vectors at time t=n;

da = a(2)-a(1);
lamM = bM*trapz(P.betaD*DH+P.betaA*AH)*da/NH; 

end