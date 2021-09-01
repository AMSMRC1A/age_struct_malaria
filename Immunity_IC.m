function [Cm,Cac,Ctot] = Immunity_IC()
global P

na = P.na;

Cm = ones(na,1); 
Cac = ones(na,1);
Ctot = P.c1*Cac+P.c2*Cm;

end