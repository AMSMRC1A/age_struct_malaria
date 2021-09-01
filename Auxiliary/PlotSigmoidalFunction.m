%% Plot Sigmoidal Probability Function

cmin = 0; % minimum value
cmax = 1; % maximum value
c = 1; % constant at half maximal
k = 2; % Hill coefficient (steepness)

x = 0:.01:5; %input

y = 1 - cmax*((1-cmin)*c^k./(c^k+x.^k)+cmin); %output

figure_setups;
plot(x,y,'linewidth',2)
hold on
plot([x(1) x(end)],[0.5 0.5],'k:')
plot([c c],[0 1],'k:')
plot([x(1) x(end)],[cmin cmin],'k:')
plot([x(1) x(end)],[cmax cmax],'k:')
hold off
text(c*1.1,0.99,'c')
text(0.9*x(end),cmin+.2*cmin,'cmin')
text(0.9*x(end),cmax+.2*cmin,'cmax')
set(gca,'fontsize',14)
xlabel('input')
ylabel('output')
