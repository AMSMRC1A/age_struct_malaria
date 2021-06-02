%% Code for Solution to Is and Phi(Is) equations
% Equations (9) and (10) in Filipe Supplemental

b = 0.25;
ds = 5;
a0 = 3;
dm = 0.25;
ks = 2;

alpha = 0:.1:60;
for EIR = [1 10 20 50 100]

    Lambda = EIR*b*(1-exp(-alpha(end)/a0));
    Is0 = 40;
    
    Is = Lambda*ds*(1+(a0*exp(-alpha/a0)-ds*exp(-alpha/ds))/(ds-a0))+Is0*exp(-alpha/dm);
    
    phi = 1./(1+(Is/Is0).^ks);


    figure(1)
    subplot(1,2,1)
    plot(alpha,phi,'linewidth',2)
    hold on
    subplot(1,2,2)
    plot(alpha,Is,'linewidth',2)
    hold on
end

subplot(1,2,1)
hold off
ylim([0 1])
title('\Phi(I_s)')
set(gca,'fontsize',16)

subplot(1,2,2)
hold off
title('I_s')
set(gca,'fontsize',16)
legend('1','10','20','50','100')
