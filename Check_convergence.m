close all
clear all
clc

age_max = 50*365; % max ages in days
dt_list = [2,5,10,20,50,60]; % solution on coarser grid
dt_r = 1; % reference solution
nrun = length(dt_list);
err = NaN(size(dt_list));

solu_r = load(['Results/solution_',num2str(dt_r),'.mat'],'P','t','Cs_t','Cs_end');
Cs_r = solu_r.Cs_end;
x_r = solu_r.P.a;

for irun = 1:nrun
    dt = dt_list(irun);
    solu = load(['Results/solution_',num2str(dt),'.mat'],'P','t','Cs_t','Cs_end');
    Cs = solu.Cs_end;
    x = solu.P.a;
    err(irun) = norm(interp1(x_r,Cs_r,x)-Cs)*solu.P.da;
end

%%
figure_setups;
loglog(dt_list,err,'*-')
polyfit(log(dt_list),log(err),1)
grid on
xlabel('dt')
ylabel('error')
title('relative order = 1.45')

figure_setups;
plot(dt_list,err,'*-')
xlabel('dt')
ylabel('error')
grid on

%%
figure_setups;
hold on
for irun = 1:nrun
    dt = dt_list(irun);
    solu = load(['Results/solution_',num2str(dt),'.mat'],'P','t','Cs_t','Cs_end');
    plot(solu.P.a,solu.Cs_end,'DisplayName',['dt=',num2str(dt)])
end
legend;
title(['~~~Immun dist at tfinal'])

figure_setups;
hold on
for irun = 1:nrun
    dt = dt_list(irun);
    solu = load(['Results/solution_',num2str(dt),'.mat'],'P','t','Cs_t','Cs_end');
    plot(solu.t,solu.Cs_t,'DisplayName',['dt=',num2str(dt)])
end
legend;
title(['~~~Total Immun in time']);