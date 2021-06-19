function X = FileRead(dir, variable, mode)
global P

fid = fopen([dir,variable,'.txt'], 'r');
switch mode
    case 'Rave' % trapz(X,1) = read and integrate in age --> time vector
        X = NaN(P.nt,1);
        for n = 1:P.nt
            X(n) = trapz(fscanf(fid, '%f ', [1,P.na]))*P.da;
        end  
    case 'Final' % X(:,end) = read only final time 
        for n = 1:P.nt-1
            fscanf(fid, '%f ', [1,P.na]);
        end  
        X = fscanf(fid, '%f ', [1,P.na]);
    case 'Entire' % X(na,nt) = read the entire matrix, large memory warning!
        X = NaN(P.na,P.nt);
        for n = 1:P.nt
            X(:,n) = fscanf(fid, '%f ', [1,P.na]);
        end
    case 'Read' % for mosquitoes only, read a time vector
        X = fscanf(fid, '%f ', [1,Inf]);
    otherwise
        error('not defined reading mode')
end
fclose(fid);
end
