function lamH = FOI_H(IM,NM,NH)
% IM NM NH(1,np)
% lamH(1,np)
global P

[bH,~] = biting_rate(NH,NM);

lamH = bH.*P.betaM.*IM./NM;

end