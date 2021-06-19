function FileSave(SH,EH,DH,AH,RH,SM,EM,IM,Cm,Cac,Cs,mode)
% mode = 'w' -- create new
% mode = 'a' -- append

fSH = fopen('Output/SH.txt', mode);
fEH = fopen('Output/EH.txt', mode);
fAH = fopen('Output/AH.txt', mode);
fDH = fopen('Output/DH.txt', mode);
fRH = fopen('Output/RH.txt', mode);
fprintf(fSH, '%12.16f ', SH); fprintf(fSH,'\n'); fclose(fSH);
fprintf(fEH, '%12.16f ', EH); fprintf(fEH,'\n'); fclose(fEH);
fprintf(fAH, '%12.16f ', AH); fprintf(fAH,'\n'); fclose(fAH);
fprintf(fDH, '%12.16f ', DH); fprintf(fDH,'\n'); fclose(fDH);
fprintf(fRH, '%12.16f ', RH); fprintf(fRH,'\n'); fclose(fRH);

fSM = fopen('Output/SM.txt', mode);
fEM = fopen('Output/EM.txt', mode);
fIM = fopen('Output/IM.txt', mode);
fprintf(fSM, '%12.16f ', SM); fprintf(fSM,'\n'); fclose(fSM);
fprintf(fEM, '%12.16f ', EM); fprintf(fEM,'\n'); fclose(fEM);
fprintf(fIM, '%12.16f ', IM); fprintf(fIM,'\n'); fclose(fIM);

fCm = fopen('Output/Cm.txt', mode);
fCac =fopen('Output/Cac.txt',mode);
fCs = fopen('Output/Cs.txt', mode);
fprintf(fCm, '%12.16f ', Cm); fprintf(fCm,'\n'); fclose(fCm);
fprintf(fCac, '%12.16f ', Cac); fprintf(fCac,'\n'); fclose(fCac);
fprintf(fCs, '%12.16f ', Cs); fprintf(fCs,'\n'); fclose(fCs);

end