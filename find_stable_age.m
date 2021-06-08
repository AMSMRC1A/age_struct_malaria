%% 
global P
global a

LP = -1;
RP = 1;
da = 1;
a_max = 50*365; % in years
a = 0:da:a_max;
%% setup the birth and mortality functions
b0 = 0;
b1 = 0.05;
b2 = 0.505;
b3 = 0.01;
b4 = 0.055;
cc = 4.6;
zz = 25;
alpha = 28;
ww = 13.5;
gH =  2*cc.*normpdf((a/365-zz)/ww).*normcdf(alpha*(a/365-zz)/ww)/ww;
P.gH = gH/365;
muH_int = (a/365)*b0 + (b1/b2)*(1-exp(-b2*a/365)) + (b3/b4)*(-1+exp(b4*a/365));
P.muH_int = muH_int/365;
%% Look for a zero in [LP,RP] to find p_hat for the stable age distribution
if stable_age_func(LP)*stable_age_func(RP)>0
    disp('No zero in that interval');
else
    p = (LP + RP)/2;
    err = abs( stable_age_func(p) );
    while err > 1e-8
        if stable_age_func(LP)*stable_age_func(p)<0
            RP = p;
        else
            LP = p;
        end
        p = (LP + RP)/2;
        err = abs( stable_age_func(p) );
    end
end
% output the final value of p_hat and the final error
disp(['Final value of p_hat = ', num2str(p)]);
disp(['Error: ', num2str(err)]);

%% Construct the stable age distribution using p_hat
Lambda = 1/(da.*trapz(exp(-p*a-P.muH_int)));
n_tilde = Lambda*exp(-p*a-P.muH_int); 

figure_setups;
plot(a,n_tilde);
xlabel('age');
ylabel('pop. density')
title('Stable Age Distribution (?)')

