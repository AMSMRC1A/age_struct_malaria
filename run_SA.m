clear all;  close all;  clc; format long;
global P lP
flag_disp = 1;

%% numerical config
tfinal = 100*365; % final time in days
age_max = 60*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);

P.a = a;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;

%% SA setting
lQ = 'R0';  % R0 RHM RMH EIR
lP_list = {'v0'}; 
% 'v0' 'bh', 'bm', 'betaM', 'betaD', 'betaA', 'muM', 'MHr', 'sigma'
Malaria_parameters_baseline;
%%
tic
if flag_disp; disp(['QOI =  ', lQ]); end
for iP = 1:length(lP_list)
    lP = lP_list{iP};   
    if flag_disp; disp(['POI =  ', lP]); end
    
    % baseline run
    if flag_disp; disp('baseline run ---'); end
    Malaria_parameters_baseline;
    P_baseline = P.(lP);
    Q_baseline = QOI_value(lQ);

    % SI index at baseline
    Qp_baseline = QOIp_value(lQ,lP);
    Qp_rescale_baseline(iP,:) = Qp_baseline'.*P_baseline'./Q_baseline;
    

    %% extended SA
    P_lower = P.([lP,'_lower']);
    P_upper = P.([lP,'_upper']);
    ngrid = 11;
    
    % allocation
    P_vals = linspace(P_lower,P_upper,ngrid)';%
    Q_vals = NaN(length(P_vals),1);
    
    for i=1:ngrid
        display(['I am working on simulation ', num2str(i)])
        P.(lP) = P_vals(i);
        Malaria_parameters_transform;
        Q_vals(i) = QOI_value(lQ);
    end
    
    %   plotting
    figure_setups; hold on
    plot(P_vals,Q_vals);
    plot(P_baseline,Q_baseline,'r*','MarkerSize',20);
    xlabel(lP)
    ylabel(lQ)
end
toc
[~,index] = sort(abs(Qp_rescale_baseline),'descend');
Qp_rescale_baseline = Qp_rescale_baseline(index);
lP_list = lP_list(index);
disp('========')
disp(lQ)
SI_index = round(Qp_rescale_baseline,2);
SS = num2str(SI_index);
for ip = 1:length(lP_list)
    display([SS(ip,:), '    ', lP_list{ip}])
end
%% generate output
% fname = sprintf('sensitivity_table_h.tex');
% Qlen = 1; Plen  = length(lP_list);
% latextable(SI_index', 'Horiz', lP_list, 'Vert', {lQ},...
%     'Hline', [0:Qlen,NaN], 'Vline', [0:Plen,NaN],...
%     'name', fname, 'format', '%.2g');
% fname = sprintf('sensitivity_table_v.tex');
% latextable(SI_index, 'Horiz', {lQ}, 'Vert', lP_list,...
%     'Hline', [0:Plen,NaN], 'Vline', [0:Qlen,NaN],...
%     'name', fname, 'format', '%.2g');