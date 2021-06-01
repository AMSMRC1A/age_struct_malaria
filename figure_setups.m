function fig = figure_setups()
% 08/31/20018 Zhuolin Qu

% Home PC single screen -- wide
% fig = figure('Position', [100, 55,1100, 600]);
% fig = figure('Position', [20, 55,1900, 500]);
% fig = figure('Position', [100, 55,1800, 900]);

% % Home PC single screen -- regular
fig = figure('Position', [100, 25,900, 700]);
% fig = figure('Position', [100, 45,900, 600]);
% Home PC single screen -- rectangle
% fig = figure('Position', [100, 55,700, 700]);

global skip_pt
global colour_r1 colour_r2 colour_r3 colour_r4
global colour_blue colour_red colour_yellow
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7

% ColorSet = varycolor(20);
% for Home PC
set(0,'defaultLineLineWidth',3,'defaultTextFontSize',30,'defaultTextFontname','CMU Serif Italic','defaultAxesFontSize',30,'defaultLineMarkerSize',8)
% for CCS imac screen
% set(0,'defaultLineLineWidth',0.5,'defaultTextFontSize',25,'defaultTextFontname','CMU Serif Italic','defaultAxesFontSize',25,'defaultLineMarkerSize',13)

% set(0,'defaultAxesLineStyleOrder',{'-','--o',':s'}) %or whatever you want
set(0,'defaultAxesColorOrder',lines(5))
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaulttextinterpreter','latex');  
set(0,'defaultLegendAutoUpdate','off');

grid on
% grid minor

box on

skip_pt = 1; % for sparser plotting

colour_r1 = [1, 0.6, 0.6]; % 80% red [255, 153, 153]
colour_r2 = [1, 0.4, 0.4]; % 70% red [255, 102, 102]
colour_r3 = [1, 0, 0]; % 50% red  red [255, 0, 0]
colour_r4 = [0.8, 0, 0]; % 40% red  red [204, 0, 0]

colour_blue = [0, 0.4470, 0.7410];
colour_red = [0.8500, 0.3250, 0.0980];
colour_yellow = [0.9290, 0.6940, 0.1250];

colour_mat1 = [0 0.4470 0.7410];
colour_mat2 = [0.8500 0.3250 0.0980];
colour_mat3 = [0.9290 0.6940 0.1250];
colour_mat4 = [0.4940 0.1840 0.5560];
colour_mat5 = [0.4660 0.6740 0.1880];
colour_mat6 = [0.3010 0.7450 0.9330];
colour_mat7 = [0.6350 0.0780 0.1840];

end




